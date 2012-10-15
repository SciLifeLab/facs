#!/usr/bin/env python

#from drass import *
import drass

print dir(drass)

drass.build("ecoli_K12.fasta", "example_filter.bloom")
drass.query("test.fastq", "example_filter.bloom")
#remove_contaminants("reference", "fastq_file")
