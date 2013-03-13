import os
import shutil
import sys
import subprocess
import unittest

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
        self._fetch_bowtie2_indices()
        
        fscreen_src = os.path.join(dirname, "fastq_screen")
        fscreen_dst = os.path.join(self.progs, "fastq_screen")
        
        if not os.path.exists(fscreen_dst):
            shutil.move(fscreen_path, self.progs)
        
        # cleanup .tar.gz & decompressed folder
        #shutil.rmtree(dirname)
        #os.remove(fname)

        cfg = open(os.path.join(self.progs, "fastq_screen.conf"), 'w')
        cfg.write(self._genconf())
        subprocess.check_call([fscreen_dst, "--conf", cfg.name])

    
    ## Aux methods for the test

    def _fetch_bowtie2_indices(self):
        genomes = []
        for ref in os.listdir(self.reference):
            # Downloads bowtie indexes genome(s)
            genomes.append(ref)
            galaxy.rsync_genomes(self.reference, genomes, ["bowtie"])

    def _genconf(self):
        for ref in os.listdir(self.reference):
            self.databases.append(("DATABASE", ref, os.path.abspath(os.path.join(ref, "indices"))))
      
        bowtie_path="bowtie2" 
        self.config = """
BOWTIE2\t\t{bowtie}
THREADS\t\t8\n
""".format(bowtie=bowtie_path)

        for db in range(len(self.databases)):
            self.config_dbs = """
{database}\t{short_name}\t{full_path}
""".format(database=self.databases[db][0], short_name=self.databases[db][1],
           full_path=self.databases[db][2])

            self.config=self.config+self.config_dbs
        return self.config

    def prepare_local_perl():
        pass
