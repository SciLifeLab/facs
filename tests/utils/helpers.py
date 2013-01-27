import os
import collections
import contextlib
import subprocess
import errno
from contextlib import contextmanager

import tempfile
from tempfile import NamedTemporaryFile
import functools
import urllib


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
            f.write(header)
            # Spike one ecoli read
            f.write(ecoli_read)

            for r in xrange(num_reads):
                # Identify reads uniquely for later debugging (task distribution, for instance)
                f.write(header + 'TASK ID: ' + str(r) + '\n')
                 
                f.write('GATTACAT' * stride + '\n')
                f.write('+' + '\n')
                f.write('arvestad' * stride + '\n')

def _download_test_files(data_dir):
    """Download required sequence and reference files.
    """

    print "Downloading reference genome files..."

    DlInfo = collections.namedtuple("ecoli", "fname dirname version")
    download_data = [DlInfo("U00096.2.fasta", "reference", None)]

    for dl in download_data:
        url = "http://togows.dbcls.jp/entry/ncbi-nucleotide/{fname}".format(fname=dl.fname)
        dirname = os.path.join(data_dir, dl.dirname)
        
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        if not os.path.exists(os.path.join(dirname, dl.fname)):
            _download_to_dir(url, dirname)

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

@contextmanager
def _make_tmp_dir():
    tmp_dir = os.environ["TMPDIR"]
    if not tmp_dir:
        home_dir = os.environ["HOME"]
        tmp_dir = os.path.join(home_dir, "tmp")
    work_dir = os.path.join(tmp_dir, "cloudbiolinux")
    if not os.path.exists(work_dir):
        subprocess.check_call(["mkdir", "-p", work_dir])
    yield work_dir
    if os.path.exists(work_dir):
        subprocess.check_call(["rm", "-rf", work_dir])

@contextmanager
def cd(path):
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)

def _fetch_and_unpack(url, need_dir=True):
    if url.startswith(("git", "svn", "hg", "cvs")):
        base = os.path.splitext(os.path.basename(url.split()[-1]))[0]
        if os.path.exists(base):
            env.safe_sudo("rm -rf {0}".format(base))
        subprocess.check_call(url)
        return base
    else:
        tar_file, dir_name, tar_cmd = _get_expected_file(url)
        print tar_file, dir_name, tar_cmd
        if not os.path.exists(tar_file):
            subprocess.check_call(["wget", "--no-check-certificate", "-O", tar_file, url])
        subprocess.check_call(["tar_cmd", "tar_file"])
        return _safe_dir_name(dir_name, need_dir)

def _get_expected_file(url):
    tar_file = os.path.split(url.split("?")[0])[-1]
    safe_tar = "--pax-option='delete=SCHILY.*,delete=LIBARCHIVE.*'"
    exts = {(".tar.gz", ".tgz") : "tar %s -xzpf" % safe_tar,
            (".tar",) : "tar %s -xpf" % safe_tar,
            (".tar.bz2",): "tar %s -xjpf" % safe_tar,
            (".zip",) : "unzip"}
    for ext_choices, tar_cmd in exts.iteritems():
        for ext in ext_choices:
            if tar_file.endswith(ext):
                return tar_file, tar_file[:-len(ext)], tar_cmd
    raise ValueError("Did not find extract command for %s" % url)

def _configure_make(env):
    subprocess.check_call(["./configure", "--disable-werror", "--prefix=", env.system_install])
    subprocess.check_call(["make", "install"])

def _get_install(url, env, make_command, post_unpack_fn=None):
    """Retrieve source from a URL and install in our system directory.
    """
    with _make_tmp_dir() as work_dir:
        with cd(work_dir):
            dir_name = _fetch_and_unpack(url)
            with cd(dir_name):
                if post_unpack_fn:
                    post_unpack_fn(env)
                make_command(env)
