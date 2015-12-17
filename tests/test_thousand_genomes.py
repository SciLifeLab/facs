import os
import sys
import errno
import glob
import facs
import unittest
import subprocess
import contextlib
import collections

from nose.plugins.attrib import attr

@attr('standard')
class ThousandGenomesTest(unittest.TestCase):
    """Build and query some simple bloom filters.
    """
    def setUp(self):
	# XXX to be factored out
        self.reference = os.path.join(os.path.dirname(__file__), "data", "reference")
        self.custom_dir = os.path.join(os.path.dirname(__file__), "data", "custom")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")
        self._install_1000g_test_files(self.custom_dir)

    def test_query_NA21137_illumina(self):
        """ Query gzip compressed fastq files
        """
        for sample in glob.glob(os.path.join(self.custom_dir, "*.fastq.gz")):
        	sample_path = os.path.join(self.custom_dir, sample)
		for ref in os.listdir(self.reference):
		    bf = os.path.join(self.bloom_dir, os.path.splitext(ref)[0]+".bloom")
		    print(sample_path, bf)
            	    facs.query(sample_path, bf)

    #XXX
    #def test_1_query_NA21137_454(self):
    #def test_1_query_NA21137_iontorrent(self):

    def _install_1000g_test_files(self, data_dir):
        """Download 1000 genomes exome data

        See ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/sequence_indices/20120522.sequence.index
        for an index of recent sequencing runs.

        Here sequencing data from individual NA21137 has (arbitrarily) been chosen for download. Sequencing was
        done at BROAD institute on a Illumina HiSeq 2000.
        """

        individual = "NA21137"
        fname = "SRR362119.filt.fastq.gz"

        base_url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/%s" % individual
        fastq_url = os.path.join(base_url, "sequence_read", fname)
        dst = os.path.join(data_dir, fname)

        if not os.path.exists(dst):
            print("downloading %s from %s" % (fname, base_url))
            cl = ["wget", fastq_url, "-O", dst]
            subprocess.check_call(cl)
