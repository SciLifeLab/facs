FACS (Fast and Accurate Classification of Sequences) C implementation
======================================================================

[![Build Status](https://travis-ci.org/SciLifeLab/facs.png?branch=master)](https://travis-ci.org/SciLifeLab/facs)
<a href="https://scan.coverity.com/projects/1599">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/1599/badge.svg"/>
</a>

WARNING: This program is under active development and this documentation might not reflect reality.
Please file a GitHub issue and we will take care of it as soon as we can.

Introduction
------------

FACS is the C implementation of a previous Perl module, please select the
<a href="https://github.com/SciLifeLab/facs/tree/perl">perl branch</a> if
you want to have a look at the old (unsupported) implementation.

Some components of this work are based in the excellent Perl `Bloom::Faster`
<a href='http://search.cpan.org/~palvaro/Bloom-Faster-1.7/lib/Bloom/Faster.pm'>implementation</a>.

Overview
--------

* 'build' is for building a bloom filter from a reference file.
It supports large genome files (>4GB), human genome, for instance.
* 'query' is for querying a fastq/fasta file against the bloom filter.
* 'remove' is for removing contamination sequences from a fastq/fasta file.


Quickstart
----------

In order to fetch the source code run:

```
$ git clone https://github.com/SciLifeLab/facs
```

For the python interface, it is highly recommended to install and run FACS under
a python virtual environment. Python virtual environments provide and isolated
environment to run your python code, solving dependency and version problems, and 
indirectly permissions. Read more about virtualenv [here](https://pypi.python.org/pypi/virtualenv).

To easily install a virtual environment you can use [virtualenv-burrito](https://github.com/brainsik/virtualenv-burrito).
Follow the instructions in the link provided in order to create a new virtual 
environment. 

Installing
----------
----------

For a standalone ```facs``` commandline tool, type: ```make```.

To compile the python bindings: ```make python```, after creating and activating the virtual environment.


Citation
--------

Henrik Stranneheim, Max Käller, Tobias Allander, Björn Andersson, Lars Arvestad, Joakim Lundeberg: Classification of DNA sequences using Bloom filters. Bioinformatics 26(13): 1595-1600 (2010)


License
-------

The code is freely available under MIT license as well as the hashing algorithm 'lookup8', which is developed by Bob Jenkins and used under MIT license.

Usage
------

Facs uses a similar commandline structure to the one found in the popular <a href="https://github.com/lh3/bwa">bwa</a>.
There are three main commands: build, query and remove. Each of them might have slightly different flags but should
behave similarly.

```
$ ./facs -h

Program: facs - Sequence analysis using bloom filters
Version: 2.0 
Contact: Enze Liu <enze.liu@scilifelab.se>

Usage:   facs <command> [options]

Command: build         build a bloom filter from a FASTA/FASTQ reference file
         query         query a bloom filter given a FASTA/FASTQ file
         remove        remove (contamination) sequences from FASTQ/FASTA file
```

For example, to build a bloom filter out of a FASTA reference genome, one should type:

```
$ ./facs build -r ecoli.fasta -o ecoli.bloom
```

That would generate a ecoli bloom filter that could be used to query a FASTQ file:

```
$ ./facs query -r ecoli.bloom -q contaminated_sample.fastq.gz -f "json"
```

Note that both plaintext fastq files and gzip-compressed files are supported transparently
to the user.

Which would return some metrics, in json format, indicating how many reads might
be contaminated with ecoli in that particular sample:

```
{
    "timestamp": "2013-03-27T11:16:21.809+0100"
    "organism": "test200.fastq"
    "bloom_filter": "eschColi_K12.bloom"
    "total_read_count": 201,
    "contaminated_reads": 1,
    "total_hits": 36,
    "contamination_rate": 0.004975,
    "p_value": 1.522929e-01
}
```

If one wishes to get `tsv` format to easily import in 
<a href="http://www.libreoffice.org/">LibreOffice.org</a> or Excel, indicate
`-f "tsv"` in the commandline, and a tsv file will be written in the local directory:

```
$ cat test200.fastq.tsv
organism    bloom_filter    total_read_count    contaminated_reads  contamination_rate
test200.fastq   eschColi_K12.bloom  201 1   0.004975
```

Finally, if one wants to remove those reads from the sample, one should run the following
command:

```
$ ./facs remove -r ecoli.bloom -q contaminated_sample.fastq
```

Output:
By using stdout and stderr, clean sequences will be stored in stdout, contaminated sequences
will be stored in stderr. They can be stored into specific files, for instance:

```
$(./facs remove -r ecoli.bloom -q contaminated_sample.fastq > clean_part.fastq ) >& contaminated_part.fastq
```

If output_path '-o' is specified, two output files will be generated:

`contaminated_sample_ecoli_contam.fastq`
`contaminated_sample_ecoli_clean.fastq`

MPI facs2.0 version
-------------------

MPI facs2.0 version can be used in multi-cpu system, for instance, a cluster, in order to take advantage 
of both multiple cores and multiple cpus at the same time.   

Usage:

First download facs package and 'make', then 'make mpi'. A unique binary file 'facs_mpi' will be generated.

```
$mpirun -np number_of_cpu ./facs_mpi -r reference_bloom_filter -q query_sequence
```
Be advised, besides openmp library, MPI facs2.0 requires MPI library (OpenMpi or Mpich, etc.)  


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
>>> facs.remove("contaminated_sample.fastq", "ecoli.bloom")
```

Update results to a database
----------------------------

FACS provides results in [JSON](http://www.json.org/) format, which eases the
storage of these results in a CouchDB instance. To do so, you need to create a
configuration file with the information for your CouchDB instance. 

The file should be named either .facsrc or .facs.cnf and should be located in 
your home directory. For system wide installations it can also be located at
/etc/facs.conf.

The format should be like this:

```
[facs]
SERVER: <your server address>
FACS_DB: <DB name>
FASTQ_SCREEN_DB: <DB name>
DECONSEQ_DB: <DB name>
USER: <username>
PASSWORD: <password>
```
