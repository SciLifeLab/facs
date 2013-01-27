#include "check.h"
#include "build.h"
#include "remove.h"
#include "remove_l.h"
#include "big_query.h"
/*------------------------------*/  
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
/*------------------------------*/ 

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1"
#endif

static int
usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: facs (Sequence decontamination using bloom filters)\n");
    fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
    fprintf(stderr, "Contact: Enze Liu <enze.liu@scilifelab.se>\n\n");
    fprintf(stderr, "Usage:   facs <command> [options]\n\n");
    fprintf(stderr, "Command: build         build a bloom filter from a FASTA reference file\n");
    fprintf(stderr, "         query         query a bloom filter given a FASTQ/FASTA file\n");
    fprintf(stderr, "         remove        remove (contamination) sequences from FASTQ/FASTA file\n");
    fprintf(stderr, "\n");
    return 1;
}


int main (int argc, char **argv) 
{
  int ret=0;
/*-------defaults-------*/ 
/*
k_mer = 21;
tole_rate = 0.8;
error_rate = 0.0005;
sampling_rate = 1;

help = 0;

prefix = NULL;
list = NULL;
ref = NULL;
source = NULL;
mode = NULL;
*/
/*-------defaults-------*/

  if (argc < 2) return usage();

  if (strcmp(argv[1], "build") == 0) ret = build_main(argc-1, argv+1);
  else if (strcmp(argv[1], "query") == 0) ret = bq_main(argc-1, argv+1);
  else if (strcmp(argv[1], "remove") == 0) ret = remove_main(argc-1, argv+1);
  else usage();
  return ret;
}
