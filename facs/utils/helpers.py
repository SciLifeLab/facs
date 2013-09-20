import os
import subprocess
import errno
import shutil
import json
import warnings
from contextlib import contextmanager

import couchdb


# Aux methods

header = '@HWI-ST188:2:1101:2751:1987#0/1'
ecoli_read = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAA"
ecoli_qual = "@Paaceeefgggfhiifghiihgiiihiiiihhhhhhhfhgcgh_fegef"

def generate_dummy_fastq(fname, num_reads, case=''):
    """ Generates simplest reads with dummy qualities
    """
    stride = 13

    if not os.path.exists(fname):
        with open(fname, "w") as f:
            # Spike one ecoli read
            f.write(header + os.linesep)
            if case == '':
                f.write(ecoli_read + os.linesep)
            elif case == '_lowercase':
                f.write(ecoli_read.lower() + os.linesep)
            elif case == '_mixedcase':
                f.write(''.join(random.choice([str.upper, str.lower])(c) for c in ecoli_read) + os.linesep)

            f.write('+' + os.linesep) # FastQ separator
            f.write(ecoli_qual + os.linesep)

            for r in xrange(num_reads):
                # Identify reads uniquely for later debugging such as
                # OpenMP parallel task distribution.
                f.write(header + 'TASK ID: ' + str(r) + os.linesep)

                f.write('GATTACAT' * stride + os.linesep)
                f.write('+' + os.linesep)
                f.write('arvestad' * stride + os.linesep)


def send_couchdb(server, db, user, passwd, doc):
    ''' Send JSON document to couchdb
    '''
    try:
        couch = couchdb.Server(server)
        couch.resource.credentials = (user, passwd)
        db = couch[db]
        db.save(json.loads(doc))
    except:
        warnings.warn("Could not connect to {server} to report test results".format(server=server))
        pass


### Software management

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

def _move_p(src, dest):
    if not os.path.exists(os.path.join(dest, os.path.split(src)[-1])):
        shutil.move(src, dest)

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
            subprocess.check_call(["rm -rf {0}".format(base)])
        subprocess.check_call(url)
        return base
    else:
        tar_file, dir_name, tar_cmd = _get_expected_file(url)
        if not os.path.exists(tar_file):
            subprocess.check_call(["wget", "--no-check-certificate", "-O", tar_file, url])
        subprocess.check_call(['tar', 'xvfz', tar_file])
        return dir_name, tar_file

def _get_expected_file(url):
    tar_file = os.path.split(url.split("?")[0])[-1]
    exts = {(".tar.gz", ".tgz"): "tar xvfz",
            (".tar",): "tar -xpf",
            (".tar.bz2",): "tar -xjpf",
            (".zip",): "unzip"}
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
