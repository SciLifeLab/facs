import os
import csv
import json
import shutil
import sys
import subprocess
import unittest

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

        self.simngs_url = 'http://www.ebi.ac.uk/goldman-srv/simNGS/current/simNGS.tgz'

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

        shutil.move(simngs, self.progs)
        shutil.move(simlib, self.progs)
        shutil.move(runfile, self.progs)

    def test_2_run_simNGS(self):
        """ Generates a synthetic library and runs with built-in simNGS runfile
        """
        reads = "100"
        simngs = os.path.join(self.progs, "simNGS")
        simlib = os.path.join(self.progs, "simLibrary")
        runfile = os.path.join(self.progs, "s_3_4x.runfile")
        ecoli = os.path.join(self.reference, "eschColi_K12", "seq", "eschColi_K12.fa")
        dst = os.path.join(self.synthetic_fastq, "simngs_{reads}.fastq".format(reads=reads))

        cl = [simlib, "-n", reads, ecoli, "|", simngs, "-o", "fastq",
              "-p", "paired", runfile, ">", dst]

        subprocess.check_call(cl)
