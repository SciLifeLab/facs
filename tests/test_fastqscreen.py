import os
import stat
import csv
import json
import glob
import shutil
import sys
import subprocess
import unittest
import datetime
from collections import defaultdict

import facs
from facs.utils import helpers, galaxy, config
from nose.plugins.attrib import attr

@attr('standard')
class FastqScreenTest(unittest.TestCase):
    """Tests against Fastq Screen, to compare performance metrics.
    """
    def setUp(self):
        self.data_dir  = os.path.join(os.path.dirname(__file__), "data")
        self.progs = os.path.join(os.path.dirname(__file__), "data", "bin")
        self.reference = os.path.join(os.path.dirname(__file__), "data", "reference")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")
        self.custom_dir = os.path.join(os.path.dirname(__file__), "data", "custom")
        self.synthetic_fastq = os.path.join(os.path.dirname(__file__), "data", "synthetic_fastq")
        self.tmp = os.path.join(os.path.dirname(__file__), "data", "tmp")

        self.fscreen_url = 'http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.4.tar.gz'

        helpers._mkdir_p(self.tmp)

        # Check if 2bit decompressor is available
        twobit_fa_path = os.path.join(self.progs, "twoBitToFa")
        if not os.path.exists(twobit_fa_path):
            galaxy.download_twoBitToFa_bin(twobit_fa_path)

        self.fastq_threads = 1
        self.databases = []
        self.results = []

    def tearDown(self):
        """ Report collated results of the tests to a remote CouchDB database.
        """
        try:
            for res in self.results:
                if config.SERVER:
                    helpers.send_couchdb(config.SERVER, config.FASTQ_SCREEN_DB, config.USERNAME, config.PASSWORD, res)
        except:
            pass

    def test_1_fetch_fastqscreen(self):
        """Downloads and installs fastq_screen locally, generates fastq_screen.conf file
        """
        # Does not work @UPPMAX, maybe some env var tweaked by the module system?
        # ... works elsewhere.
        #self.assertTrue(self._is_bowtie_present())

        dirname, fname = helpers._fetch_and_unpack(self.fscreen_url)
        self._fetch_bowtie_indices()

        fscreen_src = os.path.join(dirname, "fastq_screen")
        fscreen_dst = os.path.join(self.progs, "fastq_screen")

        if not os.path.exists(fscreen_dst):
            shutil.copy(fscreen_src, self.progs)
            shutil.copy(fscreen_src + '.conf', self.progs)

        # Install VirtualEnv Perl equivalent: cpanm
        if not os.path.exists('cpanm'):
            subprocess.check_call(['wget', 'cpanmin.us', '-O', 'cpanm'])
            os.chmod('cpanm', 0700)
            subprocess.check_call(['./cpanm', '-f', '--local-lib=~/perl5', 'local::lib'])
            subprocess.check_call(['./cpanm', '-n', '-f', 'GD::Graph::bars'])


    def test_2_run_fastq_screen(self):
        """ Runs fastq_screen tests against synthetically generated fastq files folder.
            It runs generates single threaded config files, to measure performance per-sample.
        """
        fscreen_dst = os.path.join(self.progs, "fastq_screen")
        references = glob.glob(os.path.join(self.reference, '*'))

        for fastq in glob.glob(os.path.join(self.synthetic_fastq, "*.f*q")):
            with open(os.path.join(self.progs, "fastq_screen.conf"), 'w') as cfg:
                try:
                    cfg.write(self._genconf(fastq, references.pop().split(os.path.sep)[-1], self.fastq_threads))
                except IndexError:
                    break

            fastq_path = os.path.join(self.synthetic_fastq, fastq)
            cl = ['perl', '-I', '~/perl5/lib/perl5/', '-Mlocal::lib', fscreen_dst, "--outdir", self.tmp, "--conf", cfg.name, fastq_path]
            subprocess.call(cl)

            # Process fastq_screen results format and report it in JSON
            fastq_name = os.path.basename(fastq)
            fscreen_name = os.path.splitext(fastq_name)[0]+"_screen.txt"
            fastq_screen_resfile = os.path.join(self.tmp, fscreen_name)

            if os.path.exists(fastq_screen_resfile):
                with open(fastq_screen_resfile, 'rU') as fh:
                    self.results.append(self._fastq_screen_metrics_to_json(fh, fastq_name))


    def _fastq_screen_metrics_to_json(self, in_handle, fastq_name):
        reader = csv.reader(in_handle, delimiter="\t")
        data = defaultdict(lambda: defaultdict(list))

        #Fastq_screen version: 0.4
        version = reader.next()
        # ['Library', '%Unmapped', '%One_hit_one_library', '%Multiple_hits_one_library',
        #  '%One_hit_multiple_libraries', '%Multiple_hits_multiple_libraries']
        header = reader.next()

        data['sample'] = fastq_name
        data['timestamp'] = str(datetime.datetime.utcnow())+'Z'
        data['organisms'] = []

        for row in reader:
            if not row:
                break

            organism = {}
            organism[header[0]] = row[0]
            for i in range(1,5):
                organism[header[i]] = float(row[i])

            data['organisms'].append(organism)

        return json.dumps(data)

    def _fetch_bowtie_indices(self):
        genomes = []
        for ref in os.listdir(self.reference):
            # Downloads bowtie indexes genome(s)
            genomes.append(ref)
            #XXX: parametrize for bowtie2, although it is possible that bowtie2 indices are
            # not still properly generated in the Galaxy rsync :_(
            galaxy.rsync_genomes(self.reference, genomes, ["bowtie"])

    def _genconf(self, query, reference, threads):
        bwt_index = os.path.abspath(os.path.join(self.reference, reference, "bowtie_index", reference))
        self.databases.append(("DATABASE", reference, bwt_index))

        self.config = """
BOWTIE\t\t{bowtie}
THREADS\t\t8\n
""".format(bowtie="bowtie")

        for db in range(len(self.databases)):
            self.config_dbs = """
{database}\t{short_name}\t{full_path}
""".format(database=self.databases[db][0], short_name=self.databases[db][1],
           full_path=self.databases[db][2])

            self.config=self.config+self.config_dbs
        # reset databases list, we don't want to accumulate them for
        # multithreading like before
        self.databases = []
        return self.config
