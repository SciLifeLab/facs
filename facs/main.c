#include "build.h"
#include "remove.h"
#include "query.h"
#include "info.h"

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
	   "Command: build         Build a Bloom filter from a FASTQ/FASTA reference file\n");
  fprintf (stderr,
	   "         query         Query a Bloom filter given a FASTQ/FASTA file. Large files and compressed files supported\n");
  fprintf (stderr,
	   "         remove        Remove (contamination) sequences from FASTQ/FASTA file\n");
  fprintf (stderr,
	   "         info          Print information about a Bloom filter to stdout. \n");
  fprintf (stderr, "\nUse the '-h' option to get usage information for the subcommands.\n");
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
  else if (strcmp (argv[1], "query") == 0){
    ret = bq_main (argc-1, argv+1);
  } else if (strcmp (argv[1], "remove") == 0)
    ret = remove_main (argc-1, argv+1);
  else if (strcmp (argv[1], "info") == 0)
    ret = info_main (argc-1, argv+1);
  else
    usage();

  return ret;
}
