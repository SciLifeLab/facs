import os
import re
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

        self.databases = []
        self.results = []

    def tearDown(self):
        """ Report collated results of the tests to a remote CouchDB database.
        """
        try:
            for res in self.results:
                if config.SERVER:
                    helpers.send_couchdb(config.SERVER, config.DECONSEQ_DB, config.USERNAME, config.PASSWORD, res)
        except:
            pass

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
        #""" Runs deconseq tests against synthetically generated fastq files folder.
        #    It runs generates single threaded config files, to measure performance per-sample.
        #"""
        deconseq_bin = os.path.join(self.progs, "deconseq.pl")
        references = glob.glob(os.path.join(self.reference, '*'))

        platform = os.uname()[0]

        # Use deconseq's bwa compiled in the tarball
        if platform == 'Darwin':
            bwa_bin = 'bwaMAC'
        else:
            bwa_bin = 'bwa64'
        # Exec perms not preserved in the tarball
        os.chmod(os.path.join('deconseq-standalone-0.4.3', bwa_bin), 0700)

        for fastq in glob.glob(os.path.join(self.synthetic_fastq, "*.f*q")):
            with open(os.path.join(self.progs, "DeconSeqConfig.pm"), 'w') as cfg:
                try:
                    cur_ref_path = references.pop()
                    cur_ref = cur_ref_path.split(os.path.sep)[-1]

                    cfg.write(self._genconf(cur_ref_path, cur_ref, self.tmp, self.tmp, bwa_bin))
                except IndexError:
                    break
            fastq_path = os.path.join(self.synthetic_fastq, fastq)
            cl = ['perl', deconseq_bin, "-f", fastq_path, "-dbs", cur_ref]
	    subprocess.call(cl)
	    fastq_name = os.path.basename(fastq)
            fastq_screen_resfile = os.path.join(self.tmp, fastq_name)+'.txt'
            #if os.path.exists(fastq_screen_resfile):
            with open(fastq_screen_resfile, 'w') as fh:
                 self.results.append(self._deconseq_metrics_to_json(fh, cur_ref, fastq_name))
    
    def _deconseq_metrics_to_json(self, in_handle, cur_ref, fastq_name):
        """ XXX: We should be able to find an intermediate representation/stats
            that we could use to fetch stats
        """
        data = defaultdict(lambda: defaultdict(list))
        # ['Library', '%Unmapped', '%One_hit_one_library', '%Multiple_hits_one_library',
        #  '%One_hit_multiple_libraries', '%Multiple_hits_multiple_libraries']
        #header = reader.next()
	print self._get_right_file(os.getcwd(),'cont.fq')
        data['sample'] = fastq_name
        data['timestamp'] = str(datetime.datetime.utcnow())+'Z'
	data['organisms'] = []
	organisms = {}
	contam = self._linecount_3(self._get_right_file(os.getcwd(),'cont.fq'))
	total = self._linecount_3(self._get_right_file(os.getcwd(),'clean.fq'))+contam
	organisms[cur_ref]=cur_ref
	organisms[fastq_name]=fastq_name
	organisms[total]=total
	organisms['ratio']= contam/total
	os.system("rm *cont.fq *clean.fq");
	data['organisms'].append(organisms)
	return json.dumps(data)

    def _linecount_3(self, filename):          
        count = 0
        thefile = open(filename)
        while 1:
          buffer = thefile.read(65536)
          if not buffer: break
          count += buffer.count('\n')
        return count

    def _get_right_file(self,tgt_dir,type): 
    	list_dirs = os.walk(tgt_dir) 
    	for root, dirs, files in list_dirs: 
        	for f in files: 
			if (f.find(type)>=0):
				return f
    def _fetch_bwa_indices(self):
        genomes = []
        for ref in os.listdir(self.reference):
            # Downloads bowtie indexes genome(s)
            genomes.append(ref)
            #XXX: Should only accept bwa 0.5.x, not 0.6.x indices
            galaxy.rsync_genomes(self.reference, genomes, ["bwa"])

    def _genconf(self, dbdir, ref, tmpdir, outputdir, bwa_bin):
        """Generates DeconSeq config file
        """
        dbdir = os.path.join(dbdir, "bwa_index/")
        db = os.path.basename(os.path.join(dbdir, ref+".fa"))

        self.config = """
package DeconSeqConfig;

use strict;

use constant DEBUG => 0;
use constant PRINT_STUFF => 1;
use constant VERSION => '0.4.3';
use constant VERSION_INFO => 'DeconSeq version '.VERSION;

use constant ALPHABET => 'ACGTN';

use constant DB_DIR => '%s';
use constant TMP_DIR => '%s';
use constant OUTPUT_DIR => '%s';

use constant PROG_NAME => '%s';  # should be either bwa64 or bwaMAC (based on your system architecture)
use constant PROG_DIR => 'deconseq-standalone-0.4.3/';      # should be the location of the PROG_NAME file (use './' if in the same location at the perl script)

use constant DBS => {
                      %s => {name => 'Currently testing this database in the testsuite',
                      db => '%s'}
                    };

#use constant DB_DEFAULT => '';

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
""" % (dbdir, tmpdir, outputdir, bwa_bin, ref, db)
        return self.config
