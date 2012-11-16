import os
import sys
import errno
import glob
import drass
import unittest
import subprocess
import contextlib
import collections


class DrassBasicTest(unittest.TestCase):
    """Build and query some simple bloom filters.
    """
    def setUp(self):
        self.data_dir  = os.path.join(os.path.dirname(__file__), "data")
        self.reference = os.path.join(os.path.dirname(__file__), "data", "reference")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")
        self.custom_dir = os.path.join(os.path.dirname(__file__), "data", "custom")
        self.synthetic_fastq = os.path.join(os.path.dirname(__file__), "data", "synthetic_fastq")
        self.ecoli_read = \
"""
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAA
+ 
@Paaceeefgggfhiifghiihgiiihiiiihhhhhhhfhgcgh_fegefafhhihcegbgafdbdgggceeecdd]^aWZ^Y]bba^[_b]GTXX]aOPJPSB
"""
        self.header='@HWI-ST188:2:1101:2751:1987#0/1'

        self.fastq_nreads = [1, 8, 200]

        self._mkdir_p(self.data_dir)
        self._mkdir_p(self.bloom_dir)
        self._mkdir_p(self.custom_dir)
        self._mkdir_p(self.synthetic_fastq)

        # Downloads reference genome(s)
        self._download_test_files(self.data_dir)

    def test_1_build_ref(self):
        """ Build bloom filters out of the reference genomes directory.
        """
        self.setUp()
        
        # Build bloom filter out of the reference file(s)
        for ref in os.listdir(self.reference):
            drass.build(os.path.join(self.reference, ref),
                        os.path.join(self.bloom_dir, os.path.splitext(ref)[0]+".bloom"))

    def test_2_query(self):
        """ Genereate dummy fastq files for querying purposes, query against the reference ecoli bloom filter
        """
        self.setUp()

        for nreads in self.fastq_nreads:
            test_fname = "test{}.fastq".format(nreads)
            self._generate_dummy_fastq(os.path.join(self.synthetic_fastq, test_fname), nreads)
            drass.query(os.path.join(self.synthetic_fastq, test_fname),
                        os.path.join(self.bloom_dir, "U00096.2.bloom"))
  
    def test_3_query_custom(self):
        """ Query against the uncompressed FastQ files files manually deposited in data/custom folder
        """
        self.setUp()

        for sample in glob.glob(os.path.join(self.custom_dir, "*.fastq")):
            drass.query(os.path.join(self.custom_dir, sample),
                        os.path.join(self.bloom_dir, "U00096.2.bloom"))


    def test_4_query_custom_compressed(self):
        """ Query gzip compressed fastq files
        """
        for sample in glob.glob(os.path.join(self.custom_dir, "*.fastq.gz")):
            print "Querying against compressed sample {}".format(sample)
            drass.query(os.path.join(self.custom_dir, sample),
                        os.path.join(self.bloom_dir, "U00096.2.bloom"))
  
#    def test_3_query_all_to_all_refs(self):
#        """ XXX: Query synthetic sequences against all filters?
#        """
#        pass 

# Aux methods
    def _generate_dummy_fastq(self, fname, num_reads):
        """ Generates simplest reads with dummy qualities
        """
        stride=13

        if not os.path.exists(fname):
            with open(fname, "w") as f:
                f.write(self.header)
                # Spike one ecoli read
                f.write(self.ecoli_read)

                for r in xrange(num_reads):
                    # Identify reads uniquely for later debugging (task distribution, for instance)
                    f.write(self.header + 'TASK ID: ' + str(r) + '\n')
                     
                    f.write('GATTACAT' * stride + '\n')
                    f.write('+' + '\n')
                    f.write('arvestad' * stride + '\n')

    def _download_test_files(self, data_dir):
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

    def _mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else: raise
