DASS, FACS (Fast and Accurate Classification of Sequences) C implementation
=============================================================================

Overview
--------

* 'bloom_build' is for building a bloom filter from a reference file. It Supports large genome file (>4GB),for instance, human genome.
* 'simple_check' is for checking, also supports large bloom files #XXX Checking what?
* 'simple_remove' is for removing.

Usage
------

Simply execute the scripts like './bloom_build' or './simple_remove' without any arguments will trigger an instruction page.
The whole query will be looked like these:

./bloom_build -m 1 -k 21 -e 0.005 -l reference_files_list
./simple_check -m 1 -q test.fna -l reference_files_list -t 0.8 -s 1
./simple_remove -m 1 -q test.fna -l reference_files_list -t 1

Notes
-----

* All three scripts can be executed on both Linux and Mac system. But they don't support large bloom filter building and loading on MAC system.
* If you are going to run them on any clusters in Uppmax, e.g. Kalkyl, tintin... Make sure use that when you are building a filter for a large 
genome, or when you are loading a large filter, don't directly execute them, use sbatch script instead.  Info can be found here
http://www.uppmax.uu.se/support/user-guides/kalkyl-user-guide
* These scripts support fasta and fastq format. Make sure you use the correct extension name. 
For fasta files, please use .fna or .fasta.
For fastq files, please use .fastq.
