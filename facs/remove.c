#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include<fcntl.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/time.h>
#include<sys/mman.h>
#include<sys/types.h>

#include "tool.h"
#include "bloom.h"
#include "remove.h"
#include "query.h"
#include "file_dir.h"
#ifndef __clang__
#include<omp.h>
//#include<mpi.h>
#endif

//char *clean, *contam, *clean2, *contam2;

static int
remove_usage (void)
{
  fprintf (stderr, "\nUsage: facs remove [options]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "\t-r <file>   Reference Bloom filter to query against\n");
  fprintf (stderr, "\t-q <file>   FASTA/FASTQ file containing the query reads\n");
  fprintf (stderr, "\t-t <float>  A threshold value between 0 and 1.0. Default: depends on word size (k), typically 0.4.\n");
  fprintf (stderr, "\t-o <dir>	  A directory with a '/' in the end for two output files. If a directory is not specified, clean reads and contaminated reads will be printed to stdout and stderr respectively\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "Example:\n\tfacs remove -r ecoli.bloom -q reads.fq\n");
  exit (1);
}

int remove_main (int argc, char **argv)
{
  if (argc < 2)
  {
  	return remove_usage ();
  }
/*-------defaults for bloom filter building-------*/
  int opt;
  float tole_rate = 0;
  char *ref = NULL, *list = NULL, *target_path = NULL, *source = NULL, *report_fmt = "json";
  while ((opt = getopt (argc, argv, "t:r:o:q:f:h")) != -1)
  {
      switch (opt)
      {
	case 't':
	  (optarg) && ((tole_rate = atof (optarg)), 1);
	  break;
	case 'o':
	  (optarg) && ((target_path = optarg), 1);
	  break;
	case 'q':
	  (optarg) && (source = optarg, 1);
	  break;
	case 'r':
	  (optarg) && (ref = optarg, 1);
	  break;
        case 'f': // "json", "tsv" or none
          (optarg) && (report_fmt = optarg, 1);
          break;
	case 'h':
	  return remove_usage ();
	default:
	  printf ("Unknown option: -%c\n", (char) optopt);
	  return remove_usage ();
      }
  }
  if (!target_path && !source)
  {
  	fprintf (stderr, "\nPlease, at least specify a bloom filter (-r) and a query file (-q)\n");
  	exit (-1);
  }
  char *result = query(source, ref, tole_rate, 1.000,  list, target_path, report_fmt, 'r');
  printf("%s\n",result);
  return 1;
}

/*-------------------------------------*/
void save_result (char *source, char *obj_file, char type, char *prefix, char *clean2, char *contam2)
{
  char *so = NULL, *obj = NULL;
  char *match = (char *)calloc(4 * ONE, sizeof (char)),
       *mismatch = (char *)calloc(4 * ONE, sizeof (char)),
       *so_name = (char *)calloc(4 * ONE, sizeof (char)),
       *obj_name = (char *)calloc(4 * ONE, sizeof (char));
  so = strrchr (source, '/');
  obj = strrchr (obj_file, '/');
  if (so)
  	so += 1;
  else
  	so = NULL;
  if (obj)
  	obj += 1;
  else
  	obj = NULL;
  if (so)
  	strncat (so_name, so, strrchr (source, '.') - so);
  else
  	strncat (so_name, source, strrchr (source, '.') - source);
  if (obj)
  	strncat (obj_name, obj, strrchr (obj_file, '.') - obj);
  else
  	strncat (obj_name, obj_file, strrchr (obj_file, '.') - obj_file);
  if (prefix)
  {
  	strcat (match, prefix);
  	strcat (mismatch, prefix);
  }
  else if (so)
  {
  	strncat (match, source, so - source);
  	strncat (mismatch, source, so - source);
  }
  strcat (match, so_name);
  strcat (mismatch, so_name);
  strcat (match, "_");
  strcat (mismatch, "_");
  strcat (match, obj_name);
  strcat (mismatch, obj_name);
  strcat (match, "_contam");
  strcat (mismatch, "_clean");
  if (type == '>')
  {
  	strcat (match, ".fasta");
  	strcat (mismatch, ".fasta");
  }
  else
  {
  	strcat (match, ".fastq");
  	strcat (mismatch, ".fastq");
  }
  write_result (match, contam2);
  write_result (mismatch, clean2);
  memset(contam2,0,strlen(contam2));
  memset(clean2,0,strlen(clean2));

  free (match);
  free (mismatch);
  free (so_name);
  free (obj_name);
}
/*-------------------------------------*/
