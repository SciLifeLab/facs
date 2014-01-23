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

from jinja2 import Environment, FileSystemLoader

import facs
from facs.utils import helpers, galaxy, config
from nose.plugins.attrib import attr

class DeconSeqTest(unittest.TestCase):
    """Tests against DeconSeq, to compare performance metrics.
    """
    def setUp(self):
        self.data_dir  = os.path.join(os.path.dirname(__file__), "data")
        self.progs = os.path.join(os.path.dirname(__file__), "data", "bin")
        self.reference = os.path.join(os.path.dirname(__file__), "data", "reference")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")
        self.custom_dir = os.path.join(os.path.dirname(__file__), "data", "custom")
        self.synthetic_fastq = os.path.join(os.path.dirname(__file__), "data", "synthetic_fastq")
        self.tmp = os.path.join(os.path.dirname(__file__), "data", "tmp")

        self.deconseq_url = 'http://downloads.sourceforge.net/project/deconseq/standalone/deconseq-standalone-0.4.3.tar.gz'

        helpers._mkdir_p(self.tmp)

        self.databases = []
        self.results = []

    def tearDown(self):
        """ Report collated results of the tests to a remote CouchDB database.
        """
        try:
            for res in self.results:
                if config.SERVER:
                    helpers.send_couchdb(config.SERVER, config.DECONSEQ_DB, config.USERNAME, config.PASSWORD, res, wake_up=config.WAKE)
        except:
            pass

    def test_1_fetch_deconseq(self):
        """Downloads and installs deconseq locally, generates deconseq.conf file
        """
        dirname, fname = helpers._fetch_and_unpack(self.deconseq_url)
        self._fetch_bwa_indices()

        deconseq_src = os.path.join(dirname, "deconseq.pl")
        deconseq_bin = os.path.join(self.progs, "deconseq.pl")

        if not os.path.exists(deconseq_bin):
            shutil.copy(deconseq_src, self.progs)

        # Install VirtualEnv Perl equivalent: cpanm
#        if not os.path.exists('cpanm'):
#            subprocess.check_call(['wget', 'cpanmin.us', '-O', 'cpanm'])
#            os.chmod('cpanm', 0700)
#            subprocess.check_call(['./cpanm', '-f', '--local-lib=', os.path.join(os.environ['HOME'], 'perl5'), 'local::lib'])


    def test_2_run_deconseq(self):
        """ Runs deconseq tests against synthetically generated fastq files folder.
            It runs generates single threaded config files, to measure performance per-sample.
        """
        deconseq_bin = os.path.join(self.progs, "deconseq.pl")
        references = glob.glob(os.path.join(self.reference, '*'))

        platform = os.uname()[0]

        # Use deconseq's bwa compiled in the tarball
        if platform == 'Darwin':
            bwa_bin = 'bwaMAC'
        else:
            bwa_bin = 'bwa64'

        # Exec perms not preserved in the tarball
        os.chmod(os.path.join('deconseq-standalone-0.4.3', bwa_bin), 0700)

        for fastq in glob.glob(os.path.join(self.synthetic_fastq, "*.f*q")):
            for ref in references:
                ref_shortname = os.path.basename(ref)
                with open(os.path.join(self.progs, "DeconSeqConfig.pm"), 'w') as cfg:
                    cfg.write(self._genconf(ref, ref_shortname, self.tmp, self.tmp, bwa_bin))

                fastq_path = os.path.join(self.synthetic_fastq, fastq)

                start_time = str(datetime.datetime.utcnow())+'Z'

                cl = ['perl', deconseq_bin, "-f", fastq_path, "-dbs", ref_shortname]
                subprocess.call(cl)

                end_time = str(datetime.datetime.utcnow())+'Z'

                clean = glob.glob('*_clean.fq')[0]
                contam = glob.glob('*_cont.fq')[0]

                self.results.append(self._deconseq_metrics_to_json(fastq_path, ref, clean, contam,
                                        start_time, end_time))


    def _deconseq_metrics_to_json(self, sample, ref, clean, contam, start_time, end_time):
        """ Counts contaminated and clean reads from resulting deconseq files
        """
        data = defaultdict(lambda: defaultdict(list))

        num_clean = helpers._count_lines(clean)
        num_contam = helpers._count_lines(contam)

        # FastQ format, therefore 4 lines per read
        reads_clean = int(num_clean)/4
        reads_contam = int(num_contam)/4

        if int(reads_clean) > 0:
            contamination_rate = reads_contam / reads_clean
        else: # no "clean" reads => 100% contaminated
            contamination_rate = 100

        data['contamination_rate'] = contamination_rate
        data['total_reads'] = reads_clean + reads_contam

        data['start_timestamp'] = start_time
        data['end_timestamp'] = end_time
        data['sample'] = sample
        data['reference'] = ref

        return json.dumps(data)

    def _fetch_bwa_indices(self):
        genomes = []
        for ref in os.listdir(self.reference):
            # Downloads bowtie indexes genome(s)
            genomes.append(ref)
            #XXX: Should only accept bwa 0.5.x (binary bundled within deconseq), not 0.6.x indices
            galaxy.rsync_genomes(self.reference, genomes, ["bwa"])

    def _genconf(self, dbdir, ref, tmpdir, outputdir, bwa_bin):
        """Generates DeconSeq config file
        """
        dbdir = os.path.join(dbdir, "bwa_index/")
        db = os.path.basename(os.path.join(dbdir, ref+".fa"))

        j2_env = Environment(loader=FileSystemLoader(os.path.dirname(__file__)))

        return j2_env.get_template('deconseq.j2').render(dbdir=dbdir, tmpdir=tmpdir,
                                                        outputdir=outputdir, bwa_bin=bwa_bin,
                                                        ref=ref, db=db)
