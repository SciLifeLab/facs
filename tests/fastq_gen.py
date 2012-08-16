#!/usr/bin/env python

header='@HWI-ST188:2:1101:2751:1987#0/1'
read_length=14
reads=100

ecoli_read = """
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTT 
+ 
BP\aceeefgggfhiifghiihgiiihiiiihhhhhhhfhgcgh_fegefafhhihcegbgafdbdgggceeecdd]^aWZ^Y]bba^[_b]GTXX]aOPJPS`BBBB
"""

print header
print ecoli_read

for i in xrange(reads):
    print header
    print 'GATTACAA' * read_length
    print '+'
    print 'arvestad' * read_length
