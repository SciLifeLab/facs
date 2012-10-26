#!/usr/bin/env python

#from drass import *
import drass

print dir(drass)
#drass.build("ecoli_K12.fasta", "ecoli_K12.bloom")
drass.query("test.fastq", "ecoli_K12.bloom")
#remove_contaminants("reference", "fastq_file")
