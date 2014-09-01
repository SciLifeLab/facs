#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <zlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <libgen.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>

#include "tool.h"
#include "prob.h"
#include "bloom.h"
#include "remove.h"
#include "file_dir.h"
#include "query.h"

#ifndef __clang__
  // openMP not yet ported to clang: http://www.phoronix.com/scan.php?page=news_item&px=MTI2MjU
  #include <omp.h>
#endif

char *_clean, *_contam, *_clean2, *_contam2;

static int query_usage (void)
{
  fprintf (stderr, "\nUsage: facs query [options]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "\t-r <file>    Reference Bloom filter to query against.\n");
  fprintf (stderr, "\t-q <file>    A file in FASTA/FASTQ format containing query sequences.\n");
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
  while ((opt = getopt (argc, argv, "s:t:r:o:q:f:h")) != -1)
  {
      switch (opt)
      {
	case 't':
	  tole_rate = atof(optarg);
	  break;
	case 's':
	  sampling_rate = atof(optarg); 
	  // Sampling rate is the partial proportion of a sample, or subsampling, i.e: 0.20 means take only 20% of the input file.
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
  /*
  check the format for reference bloom filter by the extension name so far, can write special tag in the filter for format checking in the future.
  */
  if (strstr(ref,".bloom") == NULL)
  {
	fprintf (stderr,"\nIncorrect bloom filter\n");
	exit(-1);
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
  int threads = 0;
  BIGCAST offset = 0;
  char *position = NULL;
  static char timestamp[40] = {0};

  // Get current timestamp, for benchmarking purposes (begin of run timestamp)
  isodate(timestamp);
  bloom *bl_2 = NEW (bloom);
  F_set *File_head = make_list (bloom_filter, list);
  /*initialization for python interface*/
  File_head->hits = 0;
  File_head->all_k = 0;
  File_head->reads_num = 0;
  File_head->reads_contam = 0;
  File_head->filename = bloom_filter;           //extra initialization for python interface
  if (load_bloom (File_head->filename, bl_2)<=0)	//load a bloom filter
	exit(-1);
  
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
  if ((zip = gzopen (query, "rb")) <= 0)
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
// XXX: Awesome will be the day when OpenMP is in OSX
#ifndef __APPLE__ 
          threads = omp_get_num_threads();
#endif
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

  /*
  mode c and r refer to contamination checking and removal function respectively. 
  The following 9 lines make sure that json/tsv output is printed after the checking 
  process, but willnot be written in stdout when running removal process.
  */
  if (target_path!=NULL || mode == 'c')
  {
  	return report(File_head, query, report_fmt, target_path, timestamp, prob_suggestion(bl_2->k_mer), threads);
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
  /*in case the file is zipped, move the satrting point to the start of real data*/
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
		pos = move_start_point (pos - 1);
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
/*move the starting point to the right position*/
char *move_start_point (char *filename)
{
  while (*filename != '\n')
  	filename--;
  filename--;			//move from \n
  while (*filename != '\n')
  	filename--;
  filename++;
  return filename;
}

void init_string(int chunk)
{
        _clean = (char *) calloc (1,chunk*sizeof (char));
        _contam = (char *) calloc (1,chunk*sizeof (char));
        _clean2 = _clean;
        _contam2 = _contam;
}

char *re_clean()
{
        return _clean2;
}

char *re_contam()
{
        return _contam2;
}

void reset_string()
{
        memset(_clean2,0,strlen(_clean2));
        memset(_contam2,0,strlen(_contam2));
        _clean = _clean2;
        _contam = _contam2;
}
/*cut the reads from the string and process them one by one*/
void read_process (bloom * bl, Queue * info, Queue * tail, F_set * File_head, float sampling_rate, float tole_rate, char mode, char fmt_type)
{
	char *start_point = info->location;
	char *next_job = NULL, *temp = NULL, *previous_point = NULL, *temp_next = NULL;
	int result = 0;
	next_job = check_fmt (info, tail, start_point, fmt_type);
	// make sure it can handle DOS and Unix format ('\r\n' and '\n')
	// XXX: what about OSX sinle '\n' ('a0' in hexa)?
	if (next_job == NULL)
		return;
	while (start_point != next_job) 
		{
		if (mode == 'c')
		{
			if (sampling_rate<1)
				temp = jump (start_point, fmt_type, sampling_rate);
			else
				temp = start_point;
		// function for fast/proportional scan
			if (start_point != temp)
			{
				start_point = temp;
				continue;
			}
		}
		// skip to the next read if needed
		#pragma omp atomic
		File_head->reads_num++;
		// atomic process for summing reads number
		previous_point = start_point;
		start_point = get_right_sp (start_point, fmt_type);
		// skip the ID line
		if (fmt_type == '@')
		{
			//identify read as fastq format read and pass it to fastq_read_check to process
			result = fastq_read_check (start_point, strchr (start_point, '\n') - start_point, 'n', bl, tole_rate, File_head);
			start_point = strchr (start_point, '\n') + 1;
			start_point = strchr (start_point, '\n') + 1;
			start_point = strchr (start_point, '\n') + 1;
		}
		else
		{
			temp_next = strchr(start_point+1,'>');
			if (temp_next == NULL)
				temp_next = next_job;
			//identify read as fasta format read and pass it to fasta_read_check to process
			result = fasta_read_check (start_point, temp_next-start_point, 'n', bl, tole_rate, File_head);
			start_point = temp_next;
		}
		if (result>0)
		{
                	 #pragma omp atomic
                         File_head->reads_contam++;
			 if (mode == 'r')
			 	{
					#pragma omp critical
					{
						memcpy(_contam,previous_point,start_point-previous_point);
						_contam+=(start_point-previous_point);
					}
				}
		}
		else
		{
			if (mode == 'r')
				{
					#pragma omp critical
					{
                                        	memcpy(_clean,previous_point,start_point-previous_point);
                                        	_clean+=(start_point-previous_point);
					}
				}
		}
	}	// outside while
}
/*generates statistic results*/
char *report(F_set *File_head, char *query, char *fmt, char *prefix, char *start_timestamp, double prob, int threads)
{
  char *abs_query_path = NULL, *abs_filter_path = NULL;
  static char buffer[800] = {0};
  static char timestamp[40] = {0};
  abs_query_path = get_abs_path(query);
  abs_filter_path = get_abs_path(File_head->filename);
  float _contamination_rate = (float) (File_head->reads_contam) / (float) (File_head->reads_num);
  double p_value = cdf(File_head->hits,get_mu(File_head->all_k,prob),get_sigma(File_head->all_k,prob));

  if(!fmt)
  {
      fprintf(stderr, "Output format not specified\n");
      exit(EXIT_FAILURE);
  } 
  else if(!strcmp(fmt, "json"))
  {
      isodate(timestamp);
      snprintf(buffer, sizeof(buffer),
"{\"begin_timestamp\": \"%s\","
"\"end_timestamp\": \"%s\","
"\"sample\": \"%s\","
"\"bloom_filter\": \"%s\","
"\"total_read_count\": %lld,"
"\"contaminated_reads\": %lld,"
"\"total_hits\": %lld,"
"\"contamination_rate\": %f,"
"\"p_value\": %e,"
"\"threads\": %d"
"}",  start_timestamp, timestamp,abs_query_path, abs_filter_path,
        File_head->reads_num, File_head->reads_contam, File_head->hits,
        _contamination_rate, p_value, threads);
  // TSV output format
  }
  else if (!strcmp(fmt, "tsv"))
  {
  	sprintf(buffer,
"sample\tbloom_filter\ttotal_read_count\t_contaminated_reads\t_contamination_rate\n"
"%s\t%s\t%lld\t%lld\t%f\t%e\n", abs_query_path , abs_filter_path,
                            File_head->reads_num, File_head->reads_contam,
                            _contamination_rate,p_value);
  }
  return buffer;
}
/*save statistic results*/
char *statistic_save (char *filename, char *prefix)
{
  char *save_file = NULL;
  int length = 0;
  if (prefix!=NULL && prefix[0]=='.')
  {
      prefix+=2;
      length = strrchr(prefix,'/')-prefix+1;
      if(length != 0 && strrchr(prefix,'/')!=NULL)
      {
      	 save_file =(char *) calloc (length, sizeof (char));
         memcpy(save_file,prefix,length);
         prefix = save_file;
         save_file = NULL;
      }
      else
      {              
         prefix = NULL;
      } 
  }
  if (prefix!=NULL)
      if (prefix[strlen(prefix)-1]=='/')
          prefix[strlen(prefix)-1]='\0'; 
  save_file = prefix_make (filename, NULL, prefix);
  if (is_dir(prefix) || prefix==NULL)
      strcat (save_file, ".info");
  if (strrchr(save_file,'/')==save_file)
      save_file++;
#ifdef DEBUG
  printf ("Basename->%s\n", filename);
  printf ("Info name->%s\n", save_file);
#endif
  return save_file;
}
/*get absolute path from a file*/
char *get_abs_path(char *filename)
{
  char *path = realpath(filename, NULL);
  if(path == NULL)
  {
        fprintf(stderr,"cannot find file with name[%s]\n", filename);
        exit(errno);
  }

  return path;
}
