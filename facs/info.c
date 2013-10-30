#include <stdio.h>
#include <stdlib.h>

#include "info.h"
#include "bloom.h"

void info_usage (void) {
  fprintf(stderr, "Usage: \n\tfacs info <file>\n");
  fprintf(stderr, "\nGet information about a FACS Bloom filter file.\n");
}


int info_main(int argc, char **argv) {
  char* bloomfile=NULL;
  bloom *bl=NULL;

  if (argc != 2) {
    info_usage();
    exit(info_bad_call);
  } else {
    bloomfile=argv[1];
    bl = NEW(bloom);
    if (bl==NULL) {
      fprintf(stderr, "Aborting: out of memory\n");
    }
  }

  load_bloom(bloomfile, bl); /* Not sure what the return vaule is. Does not seem to be an error code. /arve*/
  printf("File:\t\t%s\n", bloomfile);
  print_bloom_info(bl);

  exit(info_no_error);
}


  
