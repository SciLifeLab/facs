import os
import csv
import json
import shutil
import sys
import subprocess
import unittest
from collections import defaultdict

import facs
from facs.utils import helpers, galaxy

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

        self.databases = []

    def tearDown(self):
        # cleanup .tar.gz, decompressed folder and tmp
        fname, dirname, _ = helpers._get_expected_file(self.fscreen_url)
        try:
            shutil.rmtree(dirname)
            os.remove(fname)
            shutil.rmtree(self.tmp)
        except:
            pass

    def test_1_fetch_fastqscreen(self):
        """Downloads and installs fastq_screen locally, generates fastq_screen.conf file
        """
        #self.assertTrue(self._is_bowtie_present())

        dirname, fname = helpers._fetch_and_unpack(self.fscreen_url)
        self._fetch_bowtie_indices()

        fscreen_src = os.path.join(dirname, "fastq_screen")
        fscreen_dst = os.path.join(self.progs, "fastq_screen")

        if not os.path.exists(fscreen_dst):
            shutil.move(fscreen_path, self.progs)

        # truncates config file if present, depending on present reference genomes
        cfg = open(os.path.join(self.progs, "fastq_screen.conf"), 'w')
        cfg.write(self._genconf())
        cfg.close()

    def test_2_run_fastq_screen(self):
        """Runs fastq_screen tests against synthetically generated fastq files folder
        """
        cfg = open(os.path.join(self.progs, "fastq_screen.conf"), 'rU')
        fscreen_dst = os.path.join(self.progs, "fastq_screen")

        for fastq in os.listdir(self.synthetic_fastq):
            fastq_path = os.path.join(self.synthetic_fastq, fastq)
            cl = [fscreen_dst, "--outdir", self.tmp, "--conf", cfg.name, fastq_path]
            subprocess.call(cl)
            # Process fastq_screen results format and report it in JSON
            fastq_screen_resfile = os.path.join(self.tmp, os.path.splitext(fastq)[0]+"_screen.txt")
            if os.path.exists(fastq_screen_resfile):
                with open(fastq_screen_resfile, 'rU') as fh:
                    print self._fastq_screen_metrics_to_json(fh)

    ## Aux methods for the test
    def _is_bowtie_present(self):
        bowtie = subprocess.Popen(['which','bowtie'], shell=True, env=env,
                                  stdout=subprocess.PIPE).communicate()[0]
        # XXX: Figure out why this does behave in shell but not here
        return os.path.basename(bowtie) == "bowtie"

    def _fastq_screen_metrics_to_json(self, in_handle):
        reader = csv.reader(in_handle, delimiter="\t")
        version = reader.next()
        # ['Library', '%Unmapped', '%One_hit_one_library', '%Multiple_hits_one_library',
        #  '%One_hit_multiple_libraries', '%Multiple_hits_multiple_libraries']
        header = reader.next()
        data = defaultdict(lambda: defaultdict(dict))

        for row in reader:
            if not row:
                break
            for i in range(1,5):
                data[row[0]][header[i]] = float(row[i])
        return json.dumps(data)

    def _fetch_bowtie_indices(self):
        genomes = []
        for ref in os.listdir(self.reference):
            # Downloads bowtie indexes genome(s)
            genomes.append(ref)
            galaxy.rsync_genomes(self.reference, genomes, ["bowtie"])

    def _genconf(self):
        for ref in os.listdir(self.reference):
            bwt_index = os.path.abspath(os.path.join(self.reference, ref, "bowtie_index", ref))
            self.databases.append(("DATABASE", ref, bwt_index))

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
        return self.config

    def prepare_local_perl():
        pass
