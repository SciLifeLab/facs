FACS (Fast and Accurate Classification of Sequences) C implementation
======================================================================

[![Build Status](https://travis-ci.org/tzcoolman/DRASS.png?branch=master)](DRASS)

WARNING: This program is under active development and this documentation might not reflect reality.
Please file a GitHub issue and we will take care of it as soon as we can.

Overview
--------

* 'build' is for building a bloom filter from a reference file.
It supports large genome files (>4GB), human genome, for instance.
* 'query' is for querying a fastq/fasta file against the bloom filter.
* 'remove' is for removing contamination sequences from a fastq/fasta file.


Quickstart
----------

In order to fetch the source code, compile and run tests:

```
$ git clone https://github.com/tzcoolman/DRASS.git && cd drass && make -j8 && make tests
```

Please note that python's <a href="https://github.com/brainsik/virtualenv-burrito">virtualenv</a> is needed to run the tests.

Usage
------

Facs uses a similar commandline structure to the one found in the popular <a href="https://github.com/lh3/bwa">bwa</a>.
There are three main commands: build, query and remove. Each of them might have slightly different flags but should
behave similarly.

```
$ ./facs -h

Program: facs (Sequence decontamination using bloom filters)
Version: 0.1
Contact: Enze Liu <enze.liu@scilifelab.se>

Usage:   facs <command> [options]

Command: build         build a bloom filter from a FASTA reference file
         query         query a bloom filter given a FASTQ/FASTA file
         remove        remove (contamination) sequences from FASTQ/FASTA file
```

For example, to build a bloom filter out of a FASTA reference genome, one should type:

```
$ ./facs build -r ecoli.fasta -o ecoli.bloom
```

That would generate a ecoli bloom filter that could be used to query a FASTQ file:

```
$ ./facs query -b ecoli.bloom -r contaminated_sample.fastq.gz
```

Note that both plaintext fastq files and gzip-compressed files are supported transparently
to the user.

Which would return some metrics indicating how many reads might be contaminated with
ecoli in that particular sample:

```
{
        "total_read_count": 201,
        "contaminated_reads": 1,
        "total_hits": 90,
        "contamination_rate": 0.004975,
        "bloom_filename":"tests/data/bloom/U00096.2.bloom"
}
```

Finally, if one wants to remove those reads from the sample, one should run the following
command:

```
$ ./facs remove -b ecoli.bloom -r contaminated_sample.fastq.gz -o discarded_reads.fastq
```

Where "discarded_reads.fastq" is the reads that have been filtered out from the original
fastq file.


Python interface
----------------

A python C-Extension provides a very simple API to build, query and remove sequences,
just as described above with the plain C-based commandline.

```
$ python
Python 2.6.6 (r266:84292, Jun 18 2012, 09:57:52) 
[GCC 4.4.6 20110731 (Red Hat 4.4.6-3)] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import facs
>>> facs.build("ecoli.fasta", "ecoli.bloom")
>>> facs.query("contaminated_sample.fastq.gz", "ecoli.bloom")
>>> facs.remove("contaminated_sample.fastq.gz", "ecoli.bloom")
```


Notes
-----

* All three scripts can be executed on both Linux and Mac system. But they don't support large bloom filter building and loading on MAC system.
* FACS supports fasta and fastq formats. Make sure you use the correct extension name: .fna or .fasta and .fastq, respectively.
