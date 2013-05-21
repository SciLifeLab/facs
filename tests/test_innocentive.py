import os
import sys
import errno
import glob
import unittest
import subprocess
import contextlib
import collections

import facs

class InnocentiveTest(unittest.TestCase):
    """Test datasets from innocentive challenge @UPPMAX, not to be run outside.
       The data files are too big and do not fit Travis-CI public hosting.
    """
    def setUp(self):
        self.reference = os.path.join(os.path.dirname(__file__), "data", "reference")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")

        # Innocentive custom local paths
        self.basedir = "/proj/b2012094/Innocentive/data"
        self.testing = os.path.join(self.basedir, "Testing")
        self.example = os.path.join(self.basedir, "Example")
        self.example_results = os.path.join(self.example, "Results")

    def test_2_query_testdata(self):
        """ Query Innocentive testdata
        """
        for sample in glob.glob(os.path.join(self.testing, "*.fq")):
            print "\nQuerying against testing dataset %s" % sample
            for ref in os.listdir(self.reference):
                testset = os.path.join(self.testing, sample)
                bf = os.path.join(self.bloom_dir, os.path.splitext(ref)[0]+".bloom")
                print(testset, bf)
                facs.query(testset, bf)

    def test_3_query_exampledata(self):
        """ Query Innocentive exampledata
        """
        for sample in glob.glob(os.path.join(self.example, "*.fq")):
            print "\nQuerying against example dataset %s" % sample
            for ref in os.listdir(self.reference):
                example = os.path.join(self.example, sample)
                bf = os.path.join(self.bloom_dir, os.path.splitext(ref)[0]+".bloom")
                print(example, bf)
                facs.query(example, bf)

    def test_4_remove_host_org(self):
        """ Remove (classify) reads from hg19 into _clean and _contam files
        """
        for sample in glob.glob(os.path.join(self.example, "*.fq")):
            print "\nClassifying hg19 reads from example dataset %s" % sample
            for ref in os.listdir(self.reference):
                example = os.path.join(self.example, sample)
                bf = os.path.join(self.bloom_dir,"hg19.bloom")
                print(example, bf)
                facs.remove(example, bf)
