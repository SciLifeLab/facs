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

@attr('standard')
class FacsBasicTest(unittest.TestCase):
    """Build and query some simple bloom filters.
    """
    def setUp(self):
        self.data_dir  = os.path.join(os.path.dirname(__file__), "data")
        self.progs = os.path.join(os.path.dirname(__file__), "data", "bin")
        self.reference = os.path.join(os.path.dirname(__file__), "data", "reference")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")
        self.custom_dir = os.path.join(os.path.dirname(__file__), "data", "custom")
        self.synthetic_fastq = os.path.join(os.path.dirname(__file__), "data", "synthetic_fastq")

        self.fastq_nreads = [1, 8, 200]

        helpers._mkdir_p(self.data_dir)
        helpers._mkdir_p(self.progs)
        helpers._mkdir_p(self.reference)
        helpers._mkdir_p(self.bloom_dir)
        helpers._mkdir_p(self.custom_dir)
        helpers._mkdir_p(self.synthetic_fastq)

        # Check if 2bit decompressor is available
        twobit_fa_path = os.path.join(self.progs, "twoBitToFa")
        if not os.path.exists(twobit_fa_path):
            galaxy.download_twoBitToFa_bin(twobit_fa_path)

        # Downloads reference genome(s)
        galaxy.rsync_genomes(self.reference, ["phix", "dm3", "ecoli"], ["ucsc"], twobit_fa_path)

    def test_1_build_ref(self):
        """ Build bloom filters out of the reference genomes directory.
        """
        for ref in os.listdir(self.reference):
            org = os.path.join(self.reference, ref, "seq", ref+".fa")
            bf = os.path.join(self.bloom_dir, os.path.splitext(ref)[0]+".bloom")
            print(org, bf)
            if not os.path.exists(bf):
                facs.build(org, bf)

    def test_2_query(self):
        """ Generate dummy fastq files and query them against reference bloom filters.
        """
        for nreads in self.fastq_nreads:
            test_fname = "test%s.fastq" % nreads
            helpers.generate_dummy_fastq(os.path.join(self.synthetic_fastq,
                                                      test_fname), nreads)

        for ref in os.listdir(self.reference):
            qry = os.path.join(self.synthetic_fastq, test_fname)
            bf = os.path.join(self.bloom_dir, os.path.splitext(ref)[0]+".bloom")
            print(qry, bf)
    	    facs.query(qry, bf)

    def test_3_query_custom(self):
        """ Query against the uncompressed FastQ files files manually deposited
            in data/custom folder.
        """
        for sample in glob.glob(os.path.join(self.custom_dir, "*.fastq")):
            print "\nQuerying against uncompressed sample %s" % sample
            for ref in os.listdir(self.reference):
                facs.query(os.path.join(self.synthetic_fastq, test_fname),
                           os.path.join(self.bloom_dir, os.path.splitext(ref)[0]+".bloom"))
