#!/usr/bin/env python

#from drass import *
import drass

print dir(drass)

drass.build("../k_12.fasta", "../k_12.bloom")
#query("example_string", "example_filter.bloom", mode=1, k_mer=21, tole_rate=0.8, error_rate=0.0005, sampling_rate=1, prefix=NULL)
#remove_contaminants("reference", "fastq_file")
