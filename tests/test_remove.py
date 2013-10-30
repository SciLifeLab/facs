import os
import sys
import errno
import glob
import unittest
import subprocess
import contextlib
import collections

import facs
from facs.utils import helpers, galaxy
from nose.plugins.attrib import attr

class FacsRemoveTest(unittest.TestCase):
    """Build and query some simple bloom filters.
    """
    def setUp(self):
        self.data_dir  = os.path.join(os.path.dirname(__file__), "data")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")
        self.synthetic_fastq = os.path.join(os.path.dirname(__file__), "data", "synthetic_fastq")

    def test_2_remove(self):
        """ Remove Ecoli reads from test8 file.
        """
        qry = os.path.join(self.synthetic_fastq, "test8.fastq")
        bf = os.path.join(self.bloom_dir, "eschColi_K12.bloom")
        print(qry, bf)
        facs.remove(qry, bf)

        assert os.path.exists(os.path.join(self.synthetic_fastq, "test8_eschColi_K12_contam.fastq"))
        assert os.path.exists(os.path.join(self.synthetic_fastq, "test8_eschColi_K12_clean.fastq"))

