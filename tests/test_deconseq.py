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

class DeconSeqTest(unittest.TestCase):
    """Tests against DeconSeq, to compare performance metrics.
    """
    def setUp(self):
        self.data_dir  = os.path.join(os.path.dirname(__file__), "data")
        self.progs = os.path.join(os.path.dirname(__file__), "data", "bin")
        self.reference = os.path.join(os.path.dirname(__file__), "data", "reference")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")
        self.custom_dir = os.path.join(os.path.dirname(__file__), "data", "custom")
        self.synthetic_fastq = os.path.join(os.path.dirname(__file__), "data", "synthetic_fastq")
        self.tmp = os.path.join(os.path.dirname(__file__), "data", "tmp")

        self.deconseq_url = 'http://downloads.sourceforge.net/project/deconseq/standalone/deconseq-standalone-0.4.3.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fdeconseq%2Ffiles%2Fstandalone%2F&ts=1381593179'

        helpers._mkdir_p(self.tmp)

        # Check if 2bit decompressor is available
        twobit_fa_path = os.path.join(self.progs, "twoBitToFa")
        if not os.path.exists(twobit_fa_path):
            galaxy.download_twoBitToFa_bin(twobit_fa_path)

        self.databases = []
        self.results = []

#    def tearDown(self):
#        """ Report collated results of the tests to a remote CouchDB database.
#        """
#        try:
#            for res in self.results:
#                if config.SERVER:
#                    helpers.send_couchdb(config.SERVER, config.FASTQ_SCREEN_DB, config.USERNAME, config.PASSWORD, res)
#        except:
#            pass

    def test_1_fetch_deconseq(self):
        """Downloads and installs deconseq locally, generates deconseq.conf file
        """
        dirname, fname = helpers._fetch_and_unpack(self.deconseq_url)
        self._fetch_bwa_indices()

        deconseq_src = os.path.join(dirname, "deconseq.pl")
        deconseq_bin = os.path.join(self.progs, "deconseq.pl")

        if not os.path.exists(deconseq_bin):
            shutil.copy(deconseq_src, self.progs)

        # Install VirtualEnv Perl equivalent: cpanm
#        if not os.path.exists('cpanm'):
#            subprocess.check_call(['wget', 'cpanmin.us', '-O', 'cpanm'])
#            os.chmod('cpanm', 0700)
#            subprocess.check_call(['./cpanm', '-f', '--local-lib=', os.path.join(os.environ['HOME'], 'perl5'), 'local::lib'])


    def test_2_run_deconseq(self):
        """ Runs deconseq tests against synthetically generated fastq files folder.
            It runs generates single threaded config files, to measure performance per-sample.
        """
        deconseq_bin = os.path.join(self.progs, "deconseq")
        references = glob.glob(os.path.join(self.reference, '*'))

        for fastq in glob.glob(os.path.join(self.synthetic_fastq, "*.f*q")):
            with open(os.path.join(self.progs, "DeconSeqConfig.pm"), 'w') as cfg:
                try:
                    cfg.write(self._genconf(fastq, references.pop().split(os.path.sep)[-1], self.fastq_threads))
                except IndexError:
                    break

            fastq_path = os.path.join(self.synthetic_fastq, fastq)
            cl = ['perl', deconseq_bin, "--outdir", self.tmp, "--conf", cfg.name, fastq_path]
            subprocess.call(cl)


    def _deconseq_metrics_to_json(self, in_handle, fastq_name):
        """ XXX: We should be able to find an intermediate representation/stats
            that we could use to fetch stats
        """
        pass


    def _fetch_bwa_indices(self):
        genomes = []
        for ref in os.listdir(self.reference):
            # Downloads bowtie indexes genome(s)
            genomes.append(ref)
            #XXX: Should only accept bwa 0.5.x, not 0.6.x indices
            galaxy.rsync_genomes(self.reference, genomes, ["bwa"])

    def _genconf(self, dbdir):
        self.config = """
package DeconSeqConfig;

use strict;

use constant DEBUG => 0;
use constant PRINT_STUFF => 1;
use constant VERSION => '0.4.3';
use constant VERSION_INFO => 'DeconSeq version '.VERSION;

use constant ALPHABET => 'ACGTN';

use constant DB_DIR => '{dbdir}';
use constant TMP_DIR => '{tmpdir}';
use constant OUTPUT_DIR => '{output}';

use constant PROG_NAME => '{bwa_bin}';  # should be either bwa64 or bwaMAC (based on your system architecture)
use constant PROG_DIR => './';      # should be the location of the PROG_NAME file (use './' if in the same location at the perl script)

use constant DBS => { testdb => {name => 'Currently testing this database in the testsuite',
                             db => '{dbdir}'\}\};
use constant DB_DEFAULT => 'testdb';

#######################################################################

use base qw(Exporter);

use vars qw(@EXPORT);

@EXPORT = qw(
             DEBUG
             PRINT_STUFF
             VERSION
             VERSION_INFO
             ALPHABET
             PROG_NAME
             PROG_DIR
             DB_DIR
             TMP_DIR
             OUTPUT_DIR
             DBS
             DB_DEFAULT
             );

1;
""".format(dbdir, tmpdir, output, bwa_bin)

        return self.config
