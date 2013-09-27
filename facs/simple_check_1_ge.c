#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libgen.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>

#include "tool.h"
#include "prob.h"
#include "check.h"
#include "bloom.h"
#include "file_dir.h"

#ifndef __clang__ 
#include <omp.h>
#endif
char *clean, *contam, *clean2, *contam2;

/*save it for the possible advanced version*/

void init_string(int chunk)
{
	clean = (char *) calloc (chunk, sizeof (char));
	contam = (char *) calloc (chunk, sizeof (char));
	clean2 = clean;
	contam2 = contam;
}

char *re_clean()
{
	return clean2;
}

char *re_contam()
{
	return contam2;
}

void reset_string()
{
	memset(clean2,0,strlen(clean2));
	memset(contam2,0,strlen(contam2));
}
void read_process (bloom * bl, Queue * info, Queue * tail, F_set * File_head, float sampling_rate, float tole_rate, char mode, char fmt_type)
{
	char *start_point = info->location;
	char *next_job = NULL, *temp = NULL, *previous_point = NULL, *temp_next = NULL;
	int result = 0;
	next_job = check_fmt (info, tail, start_point, fmt_type);
	// make sure it can handle DOS and Unix format ('\r\n' and '\n')
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
						memcpy(contam,previous_point,start_point-previous_point);
						contam+=(start_point-previous_point);
					}
				}
		}
		else
		{
			if (mode == 'r')
				{
				#pragma omp critical
					{
                                        	memcpy(clean,previous_point,start_point-previous_point);
                                        	clean+=(start_point-previous_point);
					}
				}
		}
	}	// outside while
}
/*-------------------------------------*/
char *report(F_set *File_head, char *query, char *fmt, char *prefix, double prob)
{
  static char buffer[800] = {0};
  static char timestamp[40] = {0};
  float contamination_rate = (float) (File_head->reads_contam) / (float) (File_head->reads_num);
  double p_value = cdf(File_head->hits,get_mu(File_head->all_k,prob),get_sigma(File_head->all_k,prob));
  if(!fmt){
      fprintf(stderr, "Output format not specified\n");
      exit(EXIT_FAILURE);
  } else if(!strcmp(fmt, "json")) {
      isodate(timestamp);
      snprintf(buffer, sizeof(buffer),
"{\n"
"\t\"timestamp\": \"%s\",\n"
"\t\"sample\": \"%s\",\n"
"\t\"bloom_filter\": \"%s\",\n"
"\t\"total_read_count\": %lld,\n"
"\t\"contaminated_reads\": %lld,\n"
"\t\"total_hits\": %lld,\n"
"\t\"contamination_rate\": %f\n"
"\t\"p_value\": %e\n"
"}",  timestamp, query, File_head->filename,
        File_head->reads_num, File_head->reads_contam, File_head->hits,
        contamination_rate,p_value);

  // TSV output format
  } else if (!strcmp(fmt, "tsv")) {
      sprintf(buffer,
"sample\tbloom_filter\ttotal_read_count\tcontaminated_reads\tcontamination_rate\n"
"%s\t%s\t%lld\t%lld\t%f\n", query, File_head->filename,
                            File_head->reads_num, File_head->reads_contam,
                            contamination_rate,p_value);
  }
  return buffer;
}

char *statistic_save (char *filename, char *prefix)
{
  char *save_file = NULL;
  int length = 0;
  //printf ("prefix->%s\n",prefix);
  if (prefix!=NULL && prefix[0]=='.'){
      prefix+=2;
      length = strrchr(prefix,'/')-prefix+1;
      if(length != 0 && strrchr(prefix,'/')!=NULL){
         save_file =(char *) calloc (length, sizeof (char));
         memcpy(save_file,prefix,length);
         prefix = save_file;
         save_file = NULL;
      }
      else{              
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
/*-------------------------------------*/
