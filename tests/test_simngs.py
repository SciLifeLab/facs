import os
import csv
import json
import glob
import shutil
import sys
import subprocess
import unittest
from itertools import izip

import facs
from facs.utils import helpers

class SimNGSTest(unittest.TestCase):
    """ Tests against simNGS, to compare performance metrics with FACS.
    """
    def setUp(self):
        self.data_dir  = os.path.join(os.path.dirname(__file__), "data")
        self.progs = os.path.join(os.path.dirname(__file__), "data", "bin")
        self.reference = os.path.join(os.path.dirname(__file__), "data", "reference")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")
        self.custom_dir = os.path.join(os.path.dirname(__file__), "data", "custom")
        self.synthetic_fastq = os.path.join(os.path.dirname(__file__), "data", "synthetic_fastq")
        self.tmp = os.path.join(os.path.dirname(__file__), "data", "tmp")

        helpers._mkdir_p(self.data_dir)
        helpers._mkdir_p(self.progs)
        helpers._mkdir_p(self.reference)
        helpers._mkdir_p(self.bloom_dir)
        helpers._mkdir_p(self.custom_dir)
        helpers._mkdir_p(self.synthetic_fastq)
        helpers._mkdir_p(self.tmp)

        # SimNGS-specific variables
        self.simngs_url = 'http://www.ebi.ac.uk/goldman-srv/simNGS/current/simNGS.tgz'
        self.sim_reads = [100, 1000, 1000000, 10000000]

        # simNGS will generate exactly the same "random" datasets on each run
        self.sim_seed = "6666520666"
        self.simngs = os.path.join(self.progs, "simNGS")
        self.simlib = os.path.join(self.progs, "simLibrary")

        # Default Illumina error profiles for simNGS
        self.runfile = os.path.join(self.progs, "s_3_4x.runfile")

    def test_1_fetch_simNGS(self):
        """ Downloads and installs simNGS locally
        """
        dirname, fname = helpers._fetch_and_unpack(self.simngs_url)

        simngs = os.path.join(dirname, "bin", "simNGS")
        simlib = os.path.join(dirname, "bin", "simLibrary")
        runfile = os.path.join(dirname, "data", "s_3_4x.runfile")

        with helpers.cd(dirname):
            os.chdir("src")
            subprocess.check_call(["make"])

        helpers._move_p(simngs, self.progs)
        helpers._move_p(simlib, self.progs)
        helpers._move_p(runfile, self.progs)

    def test_2_run_simNGS(self):
        """ Simulates an Illumina run with simNGS read simulator
            for each organism in references directory.
        """
        # Generate N simulated reads of every organism present in "org"
        orgs = [o for o in glob.glob(os.path.join(self.reference, "*/seq/*.fa"))]

        for org in orgs:
            fa_entries = 0

            for reads in self.sim_reads:
                dst = os.path.join(self.synthetic_fastq,
                                   "simngs_{org}_{reads}.fastq".format(org=org.split(os.sep)[-3], reads=reads))

                # Synthetic file has already been generated in a previous run
                if dst is not None:
                    break

                # Determine how many FASTA "Description lines" (headers) there are
                # since simNGS will generate reads depending on that number
                with open(org, 'r') as cnt:
                    for line in cnt:
                        if '>' in line:
                            fa_entries = fa_entries+1

                with open(dst, 'w') as fh:
                    # XXX: Find a good solution for floats on reads/fa_entries
                    cl1 = [self.simlib, "--seed", self.sim_seed, "-n", str(reads/fa_entries), org]
                    cl2 = [self.simngs, "-s", self.sim_seed, "-o", "fastq", self.runfile]
                    # XXX: To be parametrized in future benchmarks (for paired end reads)
                    #cl2 = [simngs, "-o", "fastq", "-p", "paired", runfile]

                    # http://docs.python.org/2/library/subprocess.html#replacing-shell-pipeline
                    p1 = subprocess.Popen(cl1, stdout=subprocess.PIPE)
                    p2 = subprocess.Popen(cl2, stdin=p1.stdout, stdout=fh).communicate()
                    p1.stdout.close()


    def test_3_generate_mixed_dataset(self):
        """ Generates a mixed synthetic dataset of eschColi with reads[0] reads
            and dm3 with reads[1] reads.
        """
        orgs = [os.path.join(self.reference, "eschColi_K12/seq/eschColi_K12.fa"),
                os.path.join(self.reference, "dm3/seq/dm3.fa")]

        reads = [3000, 6000]

        # Will hold the real number of fastq reads after simNGS, independent of
        # organism
        lines = []

        dst = os.path.join(self.synthetic_fastq,
                           "simngs.mixed_{org1}_{org2}_{reads1}vs{reads2}.fastq".format(org1='eschColi_K12',
                                                                                        org2='dm3',
                                                                                        reads1=reads[0],
                                                                                        reads2=reads[1]))

        for org, read in izip(orgs, reads):
            with open(dst, 'a') as fh:
                cl1 = [self.simlib, "--seed", self.sim_seed, "-n", str(read), org]
                cl2 = [self.simngs, "-s", self.sim_seed, "-o", "fastq", self.runfile]

                p1 = subprocess.Popen(cl1, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(cl2, stdin=p1.stdout, stdout=fh).communicate()
                p1.stdout.close()

            lines.append(helpers._count_lines(dst))

        # Read actual FASTQ reads present in output file and rename it accordingly
        reads = [read/4 for read in lines]
        shutil.move(dst, os.path.join(self.synthetic_fastq,
                         "simngs.mixed_{org1}_{org2}_{reads1}vs{reads2}.fastq".format(org1='eschColi_K12',
                                                                                      org2='dm3',
                                                                                      reads1=reads[0],
                                                                                      reads2=reads[1]
                                                                                      )))
