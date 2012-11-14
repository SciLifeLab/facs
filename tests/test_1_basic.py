import os
#import drass
import unittest
import subprocess
import contextlib
import collections


class DrassBasicTest(unittest.TestCase):
    """Build and query some simple bloom filters.
    """
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data")
        self.bloom_dir= os.path.join(os.path.dirname(__file__), "data", "bloom")

    def test_1_build(self):
        """Test singleread fastq file.
        """
        self.setUp()
        self._install_test_files(self.data_dir)
#        print dir(drass)
#        drass.build(os.path.join(self.data_dir, "reference"))
   
   
# Aux methods
    def _install_test_files(self, data_dir):
        """Download required sequence and reference files.
        """
        DlInfo = collections.namedtuple("ecoli", "fname dirname version")
        download_data = [DlInfo("U00096.2.fasta", "reference", None)]

        for dl in download_data:
            url = "http://togows.dbcls.jp/entry/ncbi-nucleotide/{fname}".format(fname=dl.fname)
            dirname = os.path.join(data_dir, dl.dirname)
            
            if not os.path.exists(dirname):
                os.mkdir(dirname)
            if not os.path.exists(os.path.join(dirname, dl.fname)):
                self._download_to_dir(url, dirname)

    def _download_to_dir(self, url, dirname):
        print dirname
        fname = os.path.basename(url)

        cl = ["wget", url, "-O", os.path.join(dirname, fname)]
        subprocess.check_call(cl)

        # compressed tarball?
        if os.path.splitext(fname) == ".gz":
            cl = ["tar", "-xzvpf", os.path.basename(url)]

        subprocess.check_call(cl)
        os.rename(os.path.basename(dirname), dirname)
        os.remove(os.path.basename(url))
