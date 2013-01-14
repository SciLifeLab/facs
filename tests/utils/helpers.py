import os
import collections
import contextlib
import subprocess
import errno


# Aux methods

header='@HWI-ST188:2:1101:2751:1987#0/1'
ecoli_read = \
"""
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAA
+ 
@Paaceeefgggfhiifghiihgiiihiiiihhhhhhhfhgcgh_fegefafhhihcegbgafdbdgggceeecdd]^aWZ^Y]bba^[_b]GTXX]aOPJPSB
"""

def _generate_dummy_fastq(fname, num_reads):
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

def _download_test_files(data_dir):
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

def _download_to_dir(url, dirname):
    fname = os.path.basename(url)

    cl = ["wget", url, "-O", os.path.join(dirname, fname)]
    subprocess.check_call(cl)

    # compressed tarball?
    if os.path.splitext(fname) == ".gz":
        cl = ["tar", "-xzvpf", os.path.basename(url)]
        subprocess.check_call(cl)
    #os.rename(os.path.basename(dirname), dirname)
    #os.remove(os.path.basename(url))

def _mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
