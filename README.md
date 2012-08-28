DASS, FACS (Fast and Accurate Classification of Sequences) C implementation
=============================================================================

Overview
--------

* 'bloom_build' is for building a bloom filter from a reference file. It Supports large genome file (>4GB),for instance, human genome.
* 'simple_check' is for querying a fastq/fasta file against the bloom filter.
* 'simple_remove' is for removing contamination sequences from a fastq/fasta file.

Usage
------

1. ./bloom_build -r tests/data/ecoli_K12.fasta -o tests/data/ecoli.bloom
2. ./simple_check -m 1 -q tests/data/ecoli_dummy.fastq -r tests/data/ecoli.bloom

Notes
-----

* All three scripts can be executed on both Linux and Mac system. But they don't support large bloom filter building and loading on MAC system.
* If you are going to run them on any clusters in Uppmax, e.g. Kalkyl, tintin... Make sure use that when you are building a filter for a large 
genome, or when you are loading a large filter, don't directly execute them, use sbatch script instead.  Info can be found here
http://www.uppmax.uu.se/support/user-guides/kalkyl-user-guide
* These scripts support fasta and fastq format. Make sure you use the correct extension name. 
For fasta files, please use .fna or .fasta.
For fastq files, please use .fastq.
