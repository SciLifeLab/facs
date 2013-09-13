import os
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
from facs.utils import helpers, galaxy

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

        self.deconseq_url = 'http://sourceforge.net/projects/deconseq/files/standalone/deconseq-standalone-0.4.3.tar.gz/download'

        helpers._mkdir_p(self.tmp)

        # Check if 2bit decompressor is available
        twobit_fa_path = os.path.join(self.progs, "twoBitToFa")
        if not os.path.exists(twobit_fa_path):
            galaxy.download_twoBitToFa_bin(twobit_fa_path)

        self.databases = []

    def tearDown(self):
        # cleanup .tar.gz, decompressed folder and tmp
        fname, dirname, _ = helpers._get_expected_file(self.deconseq_url)
        try:
            shutil.rmtree(dirname)
            os.remove(fname)
            shutil.rmtree(self.tmp)
        except:
            pass

    def test_1_fetch_deconseq(self):
        """Downloads and installs deconseq locally, generates DeconSeq_conf.pm file
        """
        #self.assertTrue(self._is_bowtie_present())
        # Does not work @UPPMAX, maybe some env var tweaked by the module system?
        # ... works elsewhere.

        dirname, fname = helpers._fetch_and_unpack(self.deconseq_url)
        #self._fetch_bowtie_indices()

        deconseq_src = os.path.join(dirname, "deconseq")
        deconseq_dst = os.path.join(self.progs, "deconseq.pl")

        if not os.path.exists(deconseq_dst):
            shutil.move(deconseq_path, self.progs)

        # truncates config file if present, depending on present reference genomes
        cfg = open(os.path.join(self.progs, "DeconSeqConfig.pm"), 'w')
        cfg.write(self._genconf())
        cfg.close()

    def test_2_run_deconseq(self):
        """Runs deconseq tests against synthetically generated fastq files folder
        """
        cfg = open(os.path.join(self.progs, "DeconSeqConfig.pm"), 'rU')
        deconseq_dst = os.path.join(self.progs, "deconseq.pl")

        for fastq in glob.glob(os.path.join(self.synthetic_fastq, "*.f*q")):
            fastq_path = os.path.join(self.synthetic_fastq, fastq)
            cl = [deconseq_dst, "--outdir", self.tmp, "-f",fastq_path]
            subprocess.call(cl)

            # Process fastq_screen results format and report it in JSON
            fastq_name = os.path.basename(fastq)
            deconseq_name = os.path.splitext(fastq_name)[0]+"_screen.txt"
            deconseq_resfile = os.path.join(self.tmp, deconseq_name)

            if os.path.exists(deconseq_resfile):
                with open(deconseq_resfile, 'rU') as fh:
                    print self._deconseq_metrics_to_json(fh, fastq_name)

    def _deconseq_metrics_to_json(self, in_handle, fastq_name):
        reader = csv.reader(in_handle, delimiter="\t")
        data = defaultdict(lambda: defaultdict(list))

        #Fastq_screen version: 0.4
        version = reader.next()
        # ['Library', '%Unmapped', '%One_hit_one_library', '%Multiple_hits_one_library',
        #  '%One_hit_multiple_libraries', '%Multiple_hits_multiple_libraries']
        header = reader.next()

        data['sample'] = fastq_name
        data['timestamp'] = str(datetime.datetime.utcnow())+'Z'
        data['organisms'] = []

        for row in reader:
            if not row:
                break

            organism = {}
            organism[header[0]] = row[0]
            for i in range(1,5):
                organism[header[i]] = float(row[i])

            data['organisms'].append(organism)

        return json.dumps(data)

    def _fetch_bowtie_indices(self):
        genomes = []
        for ref in os.listdir(self.reference):
            # Downloads bowtie indexes genome(s)
            genomes.append(ref)
            galaxy.rsync_genomes(self.reference, genomes, ["bwa"])

    def _genconf(self):
        for ref in os.listdir(self.reference):
            bwt_index = os.path.abspath(os.path.join(self.reference, ref, "bwa_index", ref))
            self.databases.append(("DATABASE", ref, bwt_index))

        self.config = """
package DeconSeqConfig;\n
use strict;\n
use constant DEBUG => 0;\n
use constant PRINT_STUFF => 1;\n
use constant VERSION => '0.4.3';\n
use constant VERSION_INFO => 'DeconSeq version '.VERSION;\n
use constant ALPHABET => 'ACGTN';\n
use constant DB_DIR => 'db/';\n
use constant TMP_DIR => 'tmp/';\n
use constant OUTPUT_DIR => '/proj/b2012037/private/datasets/';\n
use constant PROG_NAME => 'bwa64';  use constant PROG_DIR => './';\n
use constant DBS => {\n
"""
#header
        for db in range(len(self.databases)):
            self.config_dbs = self.databases[db][0]+"=>"+"{name=>"+"\'"+"self.databases[db][0]"+"\'"+",db=>"+'\''+self.databases[db]+'\''+",parts=>1},"
#database
	self.config_end = """
};\n
use base qw(Exporter);\n
use vars qw(@EXPORT);\n
@EXPORT = qw(\n
             DEBUG\n
             PRINT_STUFF\n
             VERSION\n
             VERSION_INFO\n
             ALPHABET\n
             PROG_NAME\n
             PROG_DIR\n
             DB_DIR\n
             TMP_DIR\n
             OUTPUT_DIR\n
             DBS\n
             DB_DEFAULT\n
             );\n
1;\n
"""
            self.config=self.config+self.config_dbs+self.config_end
        return self.config

    def prepare_local_perl():
        pass
