#include "check.h"
#include "build.h"
#include "remove.h"
#include "remove_l.h"
#include "big_query.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "2.0"
#endif

static int
usage (void)
{
  fprintf (stderr, "\n");
  fprintf (stderr, "Program: facs - Sequence analysis using Bloom filters\n");
  fprintf (stderr, "Version: %s\n", PACKAGE_VERSION);
  fprintf (stderr, "Contact: Enze Liu <enze.liu@scilifelab.se>\n\n");
  fprintf (stderr, "Usage:   facs <command> [options]\n\n");
  fprintf (stderr,
	   "Command: build         build a Bloom filter from a FASTQ/FASTA reference file\n");
  fprintf (stderr,
	   "         query         query a Bloom filter given a FASTQ/FASTA file. Large files and compressed files supported\n");
  fprintf (stderr,
	   "         remove        remove (contamination) sequences from FASTQ/FASTA file\n");
  fprintf (stderr,
	   "         classify      classify reads to the most likely reference genomes\n");
  fprintf (stderr, "\n");
  return 1;
}


int
main (int argc, char **argv)
{
  int ret = 0;

  if (argc < 2)
    return usage ();

  if (strcmp (argv[1], "build") == 0)
    ret = build_main (argc-1, argv+1);
  else if (strcmp (argv[1], "query") == 0)
    ret = bq_main (argc, argv);
  else if (strcmp (argv[1], "remove") == 0)
    ret = remove_main (argc-1, argv+1);
  else if (strcmp (argv[1], "classify") == 0)
    ret = remove_l_main (argc-1, argv+1);
  else
    usage();
  return ret;
}
