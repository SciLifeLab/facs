import os
import stat
import csv
import json
import glob
import shutil
import sys
import subprocess
import unittest
import datetime
from collections import defaultdict

import facs
from facs.utils import helpers, galaxy, config
from nose.plugins.attrib import attr

try:
    from memory_profiler import memory_usage
    profile = True
except ImportError:
    profile = False

class FastqScreenTest(unittest.TestCase):
    """Tests against Fastq Screen, to compare performance metrics.
    """
    def setUp(self):
        self.data_dir  = os.path.join(os.path.dirname(__file__), "data")
        self.progs = os.path.join(os.path.dirname(__file__), "data", "bin")
        self.reference = os.path.join(os.path.dirname(__file__), "data", "reference")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")
        self.custom_dir = os.path.join(os.path.dirname(__file__), "data", "custom")
        self.synthetic_fastq = os.path.join(os.path.dirname(__file__), "data", "synthetic_fastq")
        self.tmp = os.path.join(os.path.dirname(__file__), "data", "tmp")

        self.fscreen_url = 'http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.4.2.tar.gz'

        helpers._mkdir_p(self.tmp)

        # Check if 2bit decompressor is available
        twobit_fa_path = os.path.join(self.progs, "twoBitToFa")
        if not os.path.exists(twobit_fa_path):
            galaxy.download_twoBitToFa_bin(twobit_fa_path)

        # Fastq_screen does not use OpenMP, but here we reuse the environment
        # variable from the benchmarks
        self.fastq_threads = os.environ.get('OMP_NUM_THREADS', 1)
        self.results = []

    def tearDown(self):
        """ Report collated results of the tests to a remote CouchDB database.
        """
        try:
            for res in self.results:
                if config.SERVER:
                    helpers.send_couchdb(config.SERVER, config.FASTQ_SCREEN_DB, config.USERNAME, config.PASSWORD, res, wake_up=config.WAKE)

            # remove fastq_screen files from old test runs
            self._cleanup_fastq_screen_results()
        except:
            pass

    def test_0_bowtie_versions(self):
        """ Test whether locally installed bowtie versions work, we need
            versions for provenance.
        """
        bowtie1, bowtie2 = self._bowtie_versions()

        # are the relevant version strings present in the output?
        assert "bowtie2-align" in bowtie2
        assert "bowtie" in bowtie1


    def test_1_fetch_fastqscreen(self):
        """Downloads and installs fastq_screen locally, generates fastq_screen.conf file
        """
        dirname, fname = helpers._fetch_and_unpack(self.fscreen_url)
        self._fetch_bowtie_indices()

        fscreen_src = os.path.join(dirname, "fastq_screen")
        fscreen_dst = os.path.join(self.progs, "fastq_screen")

        if not os.path.exists(fscreen_dst):
            shutil.copy(fscreen_src, self.progs)
            shutil.copy(fscreen_src + '.conf', self.progs)

        # Install VirtualEnv Perl equivalent: cpanm
        try:
            subprocess.check_call(['which', 'cpanm'])
        except:
            # Try to install fastq_screen dependencies locally, not needed in Travis
            try:
                subprocess.check_call(['wget', 'cpanmin.us', '-O', 'cpanm'])
                os.chmod('cpanm', 0700)

                perl5_local = os.path.join(os.environ['HOME'], '/perl5', 'lib', 'perl5')
                subprocess.check_call(['./cpanm', '-f', '--local-lib=', perl5_local, 'local::lib'])
            except:
                pass


    def test_2_run_fastq_screen(self):
        """ Runs fastq_screen tests against synthetically generated fastq files folder.
            It runs generates single threaded config files, to measure performance per-sample.
        """
        fscreen_dst = os.path.join(self.progs, "fastq_screen")
        references = glob.glob(os.path.join(self.reference, '*'))

        for aligner in ['bowtie', 'bowtie2']:
            print "Testing fastq_screen with {}".format(aligner)
            for fastq in glob.glob(os.path.join(self.synthetic_fastq, "*.f*q")):
                for ref in references:
                    with open(os.path.join(self.progs, "fastq_screen.conf"), 'w') as cfg:
                        cfg.write(self._genconf(fastq, ref, self.fastq_threads, aligner))

                    start_time = str(datetime.datetime.utcnow())+'Z'

                    fastq_path = os.path.join(self.synthetic_fastq, fastq)
                    cl = ['perl', '-I', os.path.join(os.environ['HOME'], "perl5/lib/perl5"), '-Mlocal::lib', fscreen_dst,
                          "--aligner", "bowtie", "--outdir", self.tmp, "--conf", cfg.name, fastq_path]
                    mem = [-1]
                    if profile:
                        mem = memory_usage(subprocess.Popen(cl), include_children=True)
                    else:
                        subprocess.call(cl)

                    end_time = str(datetime.datetime.utcnow())+'Z'

                    # Process fastq_screen results format and report it in JSON
                    fastq_name = os.path.basename(fastq)
                    fscreen_name = os.path.splitext(fastq_name)[0]+"_screen.txt"
                    fastq_screen_resfile = os.path.join(self.tmp, fscreen_name)
                    if os.path.exists(fastq_screen_resfile):
                        with open(fastq_screen_resfile, 'rU') as fh:
                            bowtie1_ver, _ = self._bowtie_versions()
                            self.results.append(self._fastq_screen_metrics_to_json(fh, \
                                    fastq_name, ref, start_time, end_time, mem, bowtie1_ver))
                        # Clean to avoid parsing the wrong results file
                        os.remove(fastq_screen_resfile)

                    self._cleanup_fastq_screen_results()


    def _fastq_screen_metrics_to_json(self, in_handle, fastq_name, ref, start_time, end_time, mem, prov):
        reader = csv.reader(in_handle, delimiter="\t")
        data = defaultdict(lambda: defaultdict(list))
        ref = os.path.basename(ref)

        #Fastq_screen version: 0.4.2
        version = reader.next()

        #['Library', '#Reads_processed', '#Unmapped', '%Unmapped',
        # '#One_hit_one_library', '%One_hit_one_library',
        # '#Multiple_hits_one_library', '%Multiple_hits_one_library',
        # '#One_hit_multiple_libraries', '%One_hit_multiple_libraries',
        # 'Multiple_hits_multiple_libraries', '%Multiple_hits_multiple_libraries']
        header = reader.next()

        data['sample'] = os.path.join(os.path.dirname(fastq_name), fastq_name)
        data['begin_timestamp'] = start_time
        data['end_timestamp'] = end_time
        data['organisms'] = []

        # Add provenance such as which version of bowtie is running
        data['bowtie_screen_version'] = prov

        for row in reader:
            # skip empty rows
            if not row:
                break

            organism = {}
            organism[header[0]] = ref
            # Go through all headers
            for i in range(1, len(header)):
                organism[header[i]] = float(row[i])
            data['organisms'].append(organism)

            # Useful to compare with other programs such as FACS or Deconseq
            print data['organisms']
            data['contamination_rate'] = data['organisms'][0]['%One_hit_one_library'] + \
                                         data['organisms'][0]['%Multiple_hits_one_library'] + \
                                         data['organisms'][0]['%One_hit_multiple_libraries']

            # Percent of mapped/unmapped should be around 100% or less
            # (XXX better way to assert this)
            assert data['contamination_rate'] + data['organisms'][0]['%Unmapped'] <= 101

            # Normalize contamination to [0, 1] values, in order to
            # make it easily comparable with those from FACS
            data['contamination_rate'] = float(data['contamination_rate']) / 100

            # reference
            data['fastq_screen_index'] = data['organisms'][0]['Library']


        # Which fastq_screen version are we running?
        data['fastq_screen_version'] = version
        # How many threads are bowtie/fastqscreen using in this test?
        data['threads'] = self.fastq_threads
        data['memory_usage'] = mem

        return json.dumps(data)

    def _fetch_bowtie_indices(self):
        genomes = []
        for ref in os.listdir(self.reference):
            # Downloads (or builds) bowtie indexes genome(s)
            genomes.append(ref)
            galaxy.rsync_genomes(self.reference, genomes, ["bowtie", "bowtie2"])


    def _bowtie_versions(self):
        """ Provenance, we should know which versions of third party tools were
            running while running the testsuite.
        """
        try:
            bowtie1_ver = subprocess.check_output(["bowtie", "--version"])
            bowtie2_ver = subprocess.check_output(["bowtie2", "--version"])
        except:
            raise RuntimeError('bowtie1 or bowtie2 seem to be missing, aborting test')

        return bowtie1_ver, bowtie2_ver

    def _genconf(self, query, reference, threads, bowtie="bowtie"):
        # The latter string (reference) shouldn't start with a slash
        bwt_index = os.path.join(self.reference, reference, bowtie+"_index", os.path.basename(reference))

        config_dbs = ""

        self.config = """
    BOWTIE\t\t{bowtie}
    THREADS\t\t{threads}\n
    """.format(bowtie=bowtie, threads=self.fastq_threads)

        config_dbs = """
    DATABASE\t{short_name}\t{full_path}\t{bowtie}
    """.format(short_name=os.path.basename(reference),
           full_path=bwt_index, bowtie=str.upper(bowtie))

        return self.config+config_dbs

    def _cleanup_fastq_screen_results(self):
        """ Delete results files across unit tests and teardown
        """
        # remove fastq_screen files from old test runs (bowtie1 vs bowtie2 have the same results files)
        shutil.rmtree(self.tmp)
        os.mkdir(self.tmp)
