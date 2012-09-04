#!/usr/bin/env python
"""
Simple dummy FastQ generator to easily test against.

Usage:
    dummy_fastq.py <num_reads> <dst_file_path>

Example:
    ./tests/fastq_dummy.py 100 ./tests/data/ecoli_dummy.fastq
"""

import os
import sys

reads=int(sys.argv[1])
dummy_fastq=os.path.join(os.path.dirname(sys.argv[2]), os.path.basename(sys.argv[2]))


ecoli_read = """
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAA
+ 
BPaaceeefgggfhiifghiihgiiihiiiihhhhhhhfhgcgh_fegefafhhihcegbgafdbdgggceeecdd]^aWZ^Y]bba^[_b]GTXX]aOPJPSB
"""
header='@HWI-ST188:2:1101:2751:1987#0/1'
stride=13

with open(dummy_fastq, "w") as f:
    f.write(header)
    f.write(ecoli_read)

    for r in xrange(reads):
        f.write(header + '\n')
        f.write('QQQQQQQQ' * stride + '\n')
        f.write('+' + '\n')
        f.write('arvestad' * stride + '\n')
