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
        
        self.bowtie_path="bowtie"
        
        helpers._mkdir_p(self.tmp)

        # Check if 2bit decompressor is available
        twobit_fa_path = os.path.join(self.progs, "twoBitToFa")
        if not os.path.exists(twobit_fa_path):
            galaxy.download_twoBitToFa_bin(twobit_fa_path)
        
        self.databases = []

    def test_1_fetch_fastqscreen(self):
        """Downloads and installs fastq_screen locally
        """
        url = 'http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.4.tar.gz'
        dirname, fname = helpers._fetch_and_unpack(url)
        self._fetch_bowtie_indices()
        
        fscreen_src = os.path.join(dirname, "fastq_screen")
        fscreen_dst = os.path.join(self.progs, "fastq_screen")
        
        if not os.path.exists(fscreen_dst):
            shutil.move(fscreen_path, self.progs)
        
        # cleanup .tar.gz & decompressed folder
        #shutil.rmtree(dirname)
        #os.remove(fname)

        cfg = open(os.path.join(self.progs, "fastq_screen.conf"), 'w')
        cfg.write(self._genconf())
        cfg.flush()

        for fastq in os.listdir(self.synthetic_fastq):
            fastq = os.path.join(self.synthetic_fastq, fastq)
            fastq_screen_resfile = os.path.join(self.tmp, os.path.splitext(fastq)[0]+"_screen.txt")
            cl = [fscreen_dst, "--outdir", self.tmp, "--conf", cfg.name, fastq]
            subprocess.call(cl)
            if os.path.exists(fastq_screen_resfile):
                print self._fastq_screen_metrics_to_json(open(fastq_screen_resfile, 'rU'))

 
    ## Aux methods for the test
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
