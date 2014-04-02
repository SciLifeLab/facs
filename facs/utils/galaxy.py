"""Retrieve indexed genomes using Galaxy's rsync server resources.

http://wiki.galaxyproject.org/Admin/Data%20Integration

Code borrowed from Brad Chapman's CloudBioLinux implementation:
    https://github.com/chapmanb/cloudbiolinux/blob/master/cloudbio/biodata/galaxy.py

Removed Fabric dependency, can only be run locally.

"""
from xml.etree import ElementTree
import glob
import os
import stat
import subprocess
from contextlib import contextmanager

# Galaxy community rsync server
server = "rsync://datacache.g2.bx.psu.edu"

# ## Compatibility definitions
index_map = {"bowtie": "bowtie_index",
             "bowtie2": "bowtie2_index",
             "bwa": "bwa_index",
             "novoalign": None,
             "ucsc": "seq",
             "seq": "sam_index"}

org_remap = {"phix": "phiX",
             "GRCh37": "hg_g1k_v37",
             "araTha_tair9": "Arabidopsis_thaliana_TAIR9",
             "araTha_tair10": "Arabidopsis_thaliana_TAIR10",
             "WS210": "ce10",
             "WS220": "ce10",
    	     "ecoli": "eschColi_K12"}

galaxy_subdirs = ["", "/microbes"]

# ## Galaxy location files

class LocCols(object):
    # Hold all possible .loc file column fields making sure the local
    # variable names match column names in Galaxy's tool_data_table_conf.xml
    def __init__(self, config, dbkey, file_path):
        self.dbkey = dbkey
        self.path = file_path
        self.value = config.get("value", dbkey)
        self.name = config.get("name", dbkey)
        self.species = config.get('species', '')
        self.index = config.get('index', 'index')
        self.formats = config.get('index', 'fastqsanger')
        self.dbkey1 = config.get('index', dbkey)
        self.dbkey2 = config.get('index', dbkey)

def _get_tool_conf(tool_name):
    """
    Parse the tool_data_table_conf.xml from installed_files subfolder and extract
    values for the 'columns' tag and 'path' parameter for the 'file' tag, returning
    those as a dict.
    """
    tool_conf = {}
    tdtc = ElementTree.parse(env.tool_data_table_conf_file)
    tables = tdtc.getiterator('table')
    for t in tables:
        if tool_name in t.attrib.get('name', ''):
            tool_conf['columns'] = t.find('columns').text.replace(' ', '').split(',')
            tool_conf['file'] = t.find('file').attrib.get('path', '')
    return tool_conf


# ## Finalize downloads

def download_twoBitToFa_bin(dst):
    platform = os.uname()[0]
    if platform == 'Darwin':
    	twobit_url = 'http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.i386/twoBitToFa'
    else:
    	twobit_url = 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa'

    subprocess.check_call(["wget", "--no-check-certificate", "{url}".format(url=twobit_url),
                           "-O", dst])
    st = os.stat(dst)
    os.chmod(dst, st.st_mode | stat.S_IEXEC)

def _finalize_index_seq(fname, binpath):
    """Convert UCSC 2bit file into fasta file.
    """
    if not binpath:
        raise ValueError("Cannot find binpath to twobittoFa")
    out_fasta = fname + ".fa"
    if not os.path.exists(out_fasta):
        subprocess.check_call([binpath, "{fname}.2bit".format(fname=fname), out_fasta])


finalize_fns = {"ucsc": _finalize_index_seq}


def _finalize_index(idx, fname, binpath):
    """Perform final processing on an rsync'ed index file if necessary.
    """
    finalize_fn = finalize_fns.get(idx)
    if finalize_fn:
        finalize_fn(fname, binpath)

# ## Retrieve data from Galaxy

def rsync_genomes(genome_dir, genomes, genome_indexes, binpath=None):
    """Top level entry point to retrieve rsync'ed indexes from Galaxy.
    """
    for gid in genomes:
        galaxy_gid = org_remap.get(gid, gid)
        indexes = _get_galaxy_genomes(galaxy_gid, genome_dir, genomes, genome_indexes)
        for idx, fname in indexes.iteritems():
            _finalize_index(idx, fname, binpath)

    return os.path.abspath(fname)

def _get_galaxy_genomes(gid, genome_dir, genomes, genome_indexes):
    """Retrieve the provided genomes and indexes from Galaxy rsync.
    """
    out = {}
    org_dir = os.path.join(genome_dir, gid)
    if not os.path.exists(org_dir):
        subprocess.check_call(['mkdir', '-p', org_dir])
    for idx in genome_indexes:
        galaxy_index_name = index_map.get(idx)
        index_file = None
        if galaxy_index_name:
            index_file = _rsync_genome_index(gid, galaxy_index_name, org_dir)
        if index_file:
            out[idx] = index_file
        else:
            print "Galaxy does not support {0} for {1}".format(idx, gid)
    return out

def _rsync_genome_index(gid, idx, org_dir):
    """Retrieve index for a genome from rsync server, returning path to files.
    """
    idx_dir = os.path.join(org_dir, idx)

    if not os.path.exists(idx_dir):
        org_rsync = None
        for subdir in galaxy_subdirs:
            test_rsync = "{server}/indexes{subdir}/{gid}/{idx}/".format(
                server=server, subdir=subdir, gid=gid, idx=idx)
            check_dir = subprocess.Popen(["rsync", "--list-only", "{server}".format(server=test_rsync)],
                                         stdout=subprocess.PIPE).communicate()[0]
            if check_dir:
                org_rsync = test_rsync
                break
        if org_rsync is None:
            # Try to build the indexes
            if not _build_index(org_dir, gid, idx.split('_')[0]):
                raise RuntimeError("Index files for {} could not be either " \
                        "downloaded or built".format(gid))

        else:
            check_dir = subprocess.Popen(["rsync", "--list-only", "{server}".format(server=org_rsync)],
                                         stdout=subprocess.PIPE).communicate()[0]

            if check_dir:
                if not os.path.exists(idx_dir):
                    subprocess.check_call(['mkdir', '-p', idx_dir])
                with cd(idx_dir):
                    subprocess.check_call(["rsync", "-avzP", org_rsync, "."])
    if os.path.exists(idx_dir):
        has_fa_ext = glob.glob("{idx_dir}/{gid}.fa*".format(idx_dir=idx_dir, gid=gid))

        ext = ".fa" if (has_fa_ext and idx not in ["seq"]) else ""
        return os.path.join(idx_dir, gid + ext)

def _build_index(genome_dir, organism, aligner):
    """Build aligner index files for the organism

    Supported aligners: BWA, Bowtie and Bowtie2.

    :param genome_dir: Directory where to store index files and where the reference genome is
    :param organism: Reference genome
    :param aligner: Aligner to build the indexes with
    :raises RuntimeError: If any error when building the indexes.
    """
    ref = os.path.join(genome_dir, 'seq', organism + '.fa')
    out_dir = os.path.join(genome_dir, aligner + '_index', organism)
    subprocess.check_call(['mkdir', '-p', out_dir])
    if aligner.find('bowtie') >= 0:
        try:
            subprocess.check_call([aligner + '-build', ref, out_dir])
            return True
        except subprocess.CalledProcessError:
            return False
    elif aligner == 'bwa':
        try:
            with cd(out_dir):
                subprocess.check_call([aligner, 'index', ref, out_dir])
                return True
        except subprocess.CalledProcessError:
            return False
    return False

@contextmanager
def cd(path):
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)
