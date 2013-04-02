import os
import sys
import glob
import subprocess
import unittest
from collections import defaultdict

import facs
from facs.utils import helpers

class FamesTest(unittest.TestCase):
    """Tests with Fames metagenomics datasets: http://fames.jgi-psf.org/Retrieve_data.html
    """
    def setUp(self):
        self.data_dir  = os.path.join(os.path.dirname(__file__), "data")
        self.reference = os.path.join(os.path.dirname(__file__), "data", "reference")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")
        self.synthetic = os.path.join(os.path.dirname(__file__), "data", "synthetic_fastq")

        self.fames_urlbase = 'ftp://ftp.jgi-psf.org/pub/JGI_data/kmavromm/fames/raw/'
        self.lowc = "simLC.seq.tar.gz"
        self.medc = "simMC.seq.tar.gz"
        self.highc = "simHC.seq.tar.gz"

    def test_1_download_low_complexity(self):
        dst = os.path.join(self.synthetic, self.lowc)

        subprocess.check_call(["wget", "{url}".format(url=self.fames_urlbase+self.lowc),
                               "-O", dst])

        with helpers.cd(self.synthetic):
            subprocess.check_call(['tar', 'xfz', self.lowc])
            os.remove(self.lowc)

    def test_2_run_low_complexity(self):
        #XXX collate the JSON output in a better way
        for sample in glob.glob(os.path.join(self.synthetic, "*.fna")):
            for ref in os.listdir(self.reference):
                facs.query(os.path.join(self.synthetic, sample),
                           os.path.join(self.bloom_dir, os.path.splitext(ref)[0]+".bloom"))

    def test_2_medium_complexity(self):
        pass

    def test_3_high_complexity(self):
        pass
