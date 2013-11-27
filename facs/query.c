#include <zlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>

#include "tool.h"
#include "prob.h"
#include "bloom.h"
#include "check.h"
#include "remove.h"
#include "file_dir.h"
#include "query.h"

#ifndef __clang__
  // openMP not yet ported to clang: http://www.phoronix.com/scan.php?page=news_item&px=MTI2MjU
  #include <omp.h>
#endif

static int query_usage (void)
{
  fprintf (stderr, "\nUsage: facs query [options]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "\t-r <file>    Reference Bloom filter to query against.\n");
  fprintf (stderr, "\t-q <file>    A file in FASTA/FASTQ format containing query sequences.\n");
  fprintf (stderr, "\t-l <file>    A file containing a list of Bloom filter files, one per line.\n");
  fprintf (stderr, "\t-t <float>   A threshold value between 0 and 1.0. Default: depends on word size (k), typically 0.4.\n");
  fprintf (stderr, "\t-f <string>  Output format for reports. Valid values are: 'json' and 'tsv'\n");
  fprintf (stderr, "\t-s <float>   Sampling rate. Setting this parameter to less than 1.0 means you only\n\t             consider a sample of reads from the query file.\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "Example:\n");
  fprintf (stderr, "\tfacs query -r hs.bloom -q reads.fq\n");
  exit(1);
}

int bq_main (int argc, char **argv)
{
  if (argc < 3)
  {
      return query_usage();
  }
/*-------defaults for bloom filter building-------*/
  int opt;
  double tole_rate = 0, sampling_rate = 1;
  char *ref = NULL, *list = NULL, *target_path = NULL, *source = NULL, *report_fmt = "json";
  // XXX: make r and l mutually exclusive
  while ((opt = getopt (argc, argv, "s:t:r:o:q:l:f:h")) != -1)
  {
      switch (opt)
      {
	case 't':
	  tole_rate = atof(optarg);
	  break;
	case 's':
	  sampling_rate = atof(optarg);
	  break;
	case 'o':
	  target_path = optarg;
	  break;
	case 'q':
	  source = optarg;
	  break;
	case 'r':
	  ref = optarg;
	  break;
	case 'l':
	  list = optarg;
	  break;
	case 'f': // "json", "tsv" or none
	  (optarg) && (report_fmt = optarg, 1);
	  break;
	case 'h':
	  return query_usage();
      	  break;
	case '?':
	  printf ("Unknown option: -%c\n", (char) optopt);
	  return query_usage();
      	  break;
      }
  }
  if (!target_path && !source)
  {
  	fprintf (stderr, "\nPlease, at least specify a bloom filter (-r) and a query file (-q)\n");
      	exit (-1);
  }
  if (target_path == NULL) 
  {
      	target_path = argv[0];
  }  //set default path, which is where the binary file is.
  char *result = query(source, ref, tole_rate, sampling_rate, list, target_path, report_fmt, 'c');
  printf("%s\n",result);
  return 1;
}

char *query (char *query, char *bloom_filter, double tole_rate, double sampling_rate, char *list, char *target_path, char *report_fmt, char mode)
{
  gzFile zip = NULL;
  char type = '@';
  int normal = 0;
  BIGCAST offset = 0;
  char *position = NULL;
  static char timestamp[40] = {0};

  // Get current timestamp, for benchmarking purposes (begin of run timestamp)
  isodate(timestamp);
  bloom *bl_2 = NEW (bloom);
  F_set *File_head = make_list (bloom_filter, list);
  /*initialize for python interface*/
  File_head->hits = 0;
  File_head->all_k = 0;
  File_head->reads_num = 0;
  File_head->reads_contam = 0;
  File_head->filename = bloom_filter;           //extra initialization for python interface
  
  load_bloom (File_head->filename, bl_2);	//load a bloom filter
  
  if (tole_rate == 0)
  {
  	tole_rate = mco_suggestion (bl_2->k_mer); // suggest an optimal match cut-off
  }
  if (mode == 'r')
  {
 	init_string(ONEG); // initialize strings for containing reads
  }
/*
  if ((get_size (query) < 2 * ONEG) && !strstr (query, ".gz") && !strstr (query, ".tar"))
        normal = 0;
  else
  {
      if ((zip = gzopen (query, "rb")) < 0)
	{
          perror ("query open error...\n");
          exit (-1);
	}
      normal = 0;
  }
*/
  if ((zip = gzopen (query, "rb")) < 0)
  {
  	fprintf(stderr, "%s\n", strerror(errno));
  	exit(EXIT_FAILURE);
  }
  if (strstr (query, ".fastq") != NULL || strstr (query, ".fq") != NULL)
  	type = '@';
  else
  	type = '>';
  if (normal == 0)
  	position = (char *) calloc (1,(ONEG+1)*sizeof (char));
  while (offset != -1)
  {
      if (normal == 1)
      {
	  position = mmaping (query);
	  offset = -1;
      }
      else
      {
	  offset = CHUNKer (zip, offset, ONEG, position, type);
      }
      Queue *head = NEW (Queue);
      head->location = NULL;
      Queue *tail = NEW (Queue);
      head->next = tail;
      Queue *head2 = head;
      get_parainfo (position, head, type);
#pragma omp parallel
      {
#pragma omp single nowait
	{
	  while (head != tail)
	    {
#pragma omp task firstprivate(head)
	      {
		if (head->location != NULL)
		{
			read_process (bl_2, head, tail, File_head, sampling_rate, tole_rate, mode, type);
		}
	      }
	      head = head->next;
	    }			// End of firstprivate
	}			// End of single - no implied barrier (nowait)
      }				// End of parallel region - implied barrier
  if (position != NULL && normal == 0)
  {
  	memset (position, 0, strlen (position));
  } 
  else if (normal == 1)
  {
	munmap (position, strlen (position));
  } 
  else
  {
 	perror ("Cannot memset, wrong position on fastq file\n");
	exit (-1);
  }
  clean_list (head2, tail);
  if (mode == 'r')
  {
	if (target_path!=NULL)
	{
      		save_result (query, File_head->filename, type, target_path, re_clean(), re_contam()); //save results into file if facs remove is called
  	}
	else
	{
		write_default(re_clean(), re_contam(), offset);
	}
	reset_string();
  }
  }				//end while
  if (normal == 0)
  {
 	bloom_destroy(bl_2);
  	gzclose(zip);
  	free (position);        //dont like file mapping, strings need to be freed in a normal way
  }
  if (target_path!=NULL)
  {
  	return report(File_head, query, report_fmt, target_path, timestamp, prob_suggestion(bl_2->k_mer));
  }
  else
  {
	char *s = "";
	return s;
  }
}

char *strrstr (char *s, char *str)
{
  char *p;
  int len = strlen (s);
  for (p = s + len - 1; p >= s; p--)
  {
  	if ((*p == *str) && (memcmp (p, str, strlen (str)) == 0))
	{
		return p;
	}
  }
  return NULL;
}

void clean_list (Queue * head, Queue * tail)
{
  Queue *element;
  while (head != tail)
  {
  	element = head->next;
  	memset (head, 0, sizeof (Queue));
  	free (head);
  	head = element;
  }
  free (tail);
}

BIGCAST CHUNKer (gzFile zip, BIGCAST offset, int chunk, char *data, char type)
{
  char c;
  char *pos = NULL;
  int length = 0;
  if (offset == 0)
  	while (offset < 10 * ONE)
    	{
  		c = gzgetc (zip);
		if (c == type)
			break;
		offset++;
  	}
  gzseek (zip, offset, SEEK_SET);
  gzread (zip, data, chunk);
  if (data != NULL)
  	length = strlen (data);
  if (length >= chunk)
  {
  	if (type == '@')
	{
		pos = strrstr (data, "\n+");
		pos = bac_2_n (pos - 1);
	}
   	else
	{
		pos = strrchr (data, '>') - 1;
	}
  }
  if (pos)
  {
  	offset += (pos - data);
  	memset (pos, 0, strlen (pos));
  }
  if (length < chunk)
  	offset = -1;
  return offset;
}

BIGCAST CHUNKgz (gzFile zip, BIGCAST offset, int chunk, char *position, char *extra, char type)
{
  memset (position, 0, chunk);
  char c, *position2 = position;
  char *x;
  int num = 0;
  if (offset == 0)
  	while (offset < 10 * ONE)
  	{
		c = gzgetc (zip);
		if (c == type)
			break;
		offset++;
  	}
  if (extra != NULL)
  {
  	memcpy (position, extra, strlen (extra));
  	position += strlen (extra);
  }
  free (extra);
  while (((c = gzgetc (zip)) != EOF) && (num < chunk))
  {
  	*position = c;
  	position++;
  	num++;
  }
  x = strrstr (position2, "\n@");
  extra=(char *)calloc(1,(position-x+1)*sizeof(char));
  memcpy(x,extra,position-x+1);
  offset+=(position-x+1);
  return offset;
}

char *bac_2_n (char *filename)
{
  while (*filename != '\n')
  	filename--;
  filename--;			//move from \n
  while (*filename != '\n')
  	filename--;
  filename++;
  return filename;
}

