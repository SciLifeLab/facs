#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>
#include <sys/time.h>

#ifndef __clang__
#include <omp.h>
#endif

#include "tool.h"
#include "bloom.h"
#include "file_dir.h"
/*
In the following functions, there are two different read formats and two read scanning mode
'q' means fastq format; 'a' means fasta format; 'n' means normal mode; 'r' means reverse compliment mode.
One read will be scaned in normal mode first and reverse complement mode next. Only if the read has
been identified as a hit read, the reverse complement mode will be skipped.
*/

void isodate(char* buf) {
    /* Borrowed from: https://raw.github.com/jordansissel/experiments/bd58235b99f608472212a5933b52fca9cf1cac8d/c/time/iso8601.c */
    struct timeval tv;
    struct tm tm;
    char timestamp[] = "YYYY-MM-ddTHH:mm:ss.SSS+0000";

    /* Get the current time at high precision; could also use clock_gettime() for
     * even higher precision times if we want it. */

    gettimeofday(&tv, NULL);

    /* convert to time to 'struct tm' for use with strftime */

    localtime_r(&tv.tv_sec, &tm);

    /* format the time */

    strftime(timestamp, sizeof(timestamp), "%Y-%m-%dT%H:%M:%S.000%z", &tm);

    /* but, since strftime() can't subsecond precision, we have to hack it
     * in manually. '20' is the string offset of the subsecond value in our
     * timestamp string. Also, because sprintf always writes a null, we have to
     * write the subsecond value as well as the rest of the string already there.
     */

    sprintf(timestamp + 20, "%03d%s", tv.tv_usec / 1000, timestamp + 23);
    sprintf(buf, "%s", timestamp);
}
/*sub function for quick pass*/
int total_subscan (bloom *bl, F_set *File_head, char *begin, char *start_point, int read_length, int true_length, float tole_rate, char mode, char type)
{
	
	if (mode == 'n')
	{
		#pragma omp atomic
		File_head->all_k += (true_length);
	}
	
	int result = 0;
	while (read_length > 0)
        {
                if (read_length >= bl->k_mer)
                {
                        read_length -= bl->k_mer;
                }
                else
                {
                        start_point -= (bl->k_mer-read_length);
                        read_length = 0;
                }
                if (bloom_check (bl, start_point))
                {
                        result = total_full_check (bl, begin, true_length, tole_rate, File_head);
                        if (result > 0)
                        {
                                if (mode == 'r') 
                                        free(begin);
                                return result;
                        }
                        else if (mode == 'n')
                                break;
                }
                start_point+=bl->k_mer;
        }               
        if (mode == 'r')
        {
                free(begin); 
                return 0;
        }
        else
        {
		if (type=='q')
                	return fastq_read_check (begin, true_length, 'r', bl, tole_rate, File_head);
        	else
			return fasta_read_check (begin, true_length, 'r', bl, tole_rate, File_head);
	}
}
/*quick pass for fastq reads using k-mer and 0 overlap
if one hit exists, pass the read to total_subscan to do full check and return the value, other wise return 0
*/
int fastq_read_check (char *begin, int length, char mode, bloom * bl, float tole_rate, F_set * File_head)
{
/*
begin is the start point of the read in the string. Hence length is read length. mode contains 'r' and 'n'
as mentioned above. bl is a structor containing bloom filter. tole_rate is the threshold/match cutoff  value 
File_head is a structor containing all the info for the query, including statistic info.
Same principle applies to fasta_read_check function.
*/
	if (mode == 'r') // make a copy of the read for reverse compliment process
	{
		char *re_compliment = (char *) calloc (length+1, sizeof (char));
		re_compliment[length]='\0';
		memcpy(re_compliment, begin, length);
		begin = re_compliment;
		rev_trans (begin,length);
	}
	// initialization
	int read_length = length;
	char *start_point = begin;
	if (mode == 'n')
	{
		normal_lower(start_point,length); //normalize the whole read tddo the lower case
	}
	return total_subscan (bl, File_head, begin, start_point, read_length, length, tole_rate, mode, 'q');
}
/*full check for fastq or fasta sequence with k-mer and k-1 overlap
return a value that larger than match cut off, otherwise, return 0
*/
int total_full_check (bloom * bl, char *start_point, int length, float tole_rate, F_set * File_head)
{
	int read_length = length, count = 0, match_s = 0, mark = 1, prev = 0, conse = 0, match_time = 0;
	float result;
	while (read_length >= bl->k_mer)
	{
      		if (count >= bl->k_mer)
		{
	  		mark = 1;
	  		count = 0;
		}
		if (bloom_check (bl, start_point))
		{	
			match_time++;
			if (prev == 1)
				conse++;
			else
			{
				conse += bl->k_mer;
				prev = 1;
			}
			if (mark == 1)
			{
				match_s += (bl->k_mer - 1);
				mark = 0;
			}
			else
			match_s++;
		}
		else
		{
			prev = 0;
		}
			count++;
		start_point++;
		read_length--;
	}				// end while
	result = (float)(match_time*bl->k_mer+conse)/(float)(bl->k_mer*length+length-bl->dx);
	//result = (float) (match_time * bl->k_mer + conse) / (float) (length * bl->k_mer - 2 * bl->dx + length - bl->k_mer + 1);
	#pragma omp atomic
	File_head->hits += match_time;
	if (result >= tole_rate)
		return match_s;
	else
		return 0;
}
/*fasta read quick check using k-mer and 0 overlap*/
int fasta_read_check (char *begin, int length, char mode, bloom * bl, float tole_rate, F_set * File_head)
{
	// skip id line
	char *start_point = NULL;
	//int result; 
	int true_length = 0, read_length = 0;
	if (!begin || *begin == '>')
  		return 1;
	// in case the read is empty
	if (mode == 'n')
		start_point = fa_count (begin, length);
	else
		start_point = begin;
	true_length = strlen(start_point);
	read_length = true_length;
       	if (mode == 'r') // make a copy of the read for reverse compliment process
       	{
                rev_trans (start_point,true_length); // reverse compliment process
        }
	if (mode == 'n')
		begin = start_point;
        normal_lower(start_point,true_length); 
	return total_subscan (bl, File_head, begin, start_point, read_length, true_length, tole_rate, mode, 'a');
}
/*Parallel job distribution*/
int get_parainfo (char *full, Queue * head, char type)
{
#ifdef DEBUG
  printf ("distributing...\n");
#endif
      char *previous = NULL;
	  char *temp = full;
#ifndef __clang__
	  int cores = omp_get_num_procs ();
#else
	  int cores = 1;
#endif
      short add = 0;
      int offset = 0;
      Queue *pos = head;
      int length = 0;
      if (full != NULL)
      {
          offset = strlen(full) / cores;
      }
      if (type == '>')
      {
              for (add = 0; add < cores; add++)
              {
                  Queue *x = NEW (Queue);
		  
                  if (add == 0 && *full != '>')
			{
                  	temp = strchr (full, '>');	//drop the possible fragment
			}                
		  if (add != 0)
                  	temp = strchr (full + offset * add, '>');
                  x->location = temp;
                  x->number = &add;
                  x->next = pos->next;
                  pos->next = x;
                  pos = pos->next;
              }
      } 
      else
      {
	      for (add = 0; add < cores; add++) 
              {
              	char *tx = strchr(full,'\n');
              	length = strchr(tx+1,'\n')-(tx+1);
	      	Queue *x = NEW (Queue);
              	x->location = NULL;
		if (add != 0)
			{
                	temp = fastq_relocate(full, offset*add, length);   
			}    
              	if (previous!=temp)
		{
                	previous = temp;
                	x->location = temp;
                	x->number = &add;
                	x->next = pos->next;
                	pos->next = x;
                	pos = pos->next;
                }
	      }
      }
  return type;
}

/*reads skipping process for proportional check*/
char *jump (char *target, char type, float sampling_rate)
{
  float seed = rand () % 10;

  if (seed >= (float) sampling_rate * 10)
    {

      char *point;

      if (type == '>')
	point = strchr (target + 1, '>');	//point to >
      else
	{
	  //point = strstr (target + 1, "\n+") + 1;	//point to +
          //for (x=0;x<4;x++)
        point = strchr (target, '\n') + 1; 
        point = strchr (point, '\n') + 1; 
        point = strchr (point, '\n') + 1; 
	point = strchr (point, '\n') + 1;	//point to quality line
	}
      if (point)
	target = point;
    }
  return target;
}
/*relocate the starting points (correct @ positions) for fastq files*/
char *fastq_relocate (char *data, int offset, int length)
{
  char *target = NULL;
  int current_length = 0, read_length = 0;
  if (data != NULL && offset != 0)
    {
      target = strstr (data + offset, "\n+");
      if (!target)
	return NULL;
      else
	{
	  current_length = strchr (target + 1, '\n') - target + 1;
	  read_length = fq_read_length (target - 1);
	  if (read_length != current_length)
	    target = strchr (target + 1, '\n') + 1;
	  if (target != NULL)
	    target = strchr (target + 1, '\n') + 1;
	}
    }
  return target;
}


/*get read length for fastq file*/
int fq_read_length (char *data)
{
	char *origin = data;
	while (*data != '\n')
		data--;
	return origin - data;
}
/*check the head of the file and see if it is standard*/
char *check_fmt (Queue *info, Queue *tail, char *start_point, char type)
{
	char *next_job = NULL;
        if(info->location[0] != type)
        {
                return next_job;
        }
        else if(info->next != tail && info->next->location != NULL)
        {
                next_job = info->next->location;
        }
        else
        { // XXX: Does this account for OSX-style newlines? next_job is always NULL on OSX fastq files.
                next_job = strchr (start_point, '\0');
                if (next_job[-1] == '\n' && next_job[-2] == '\n')
                        next_job -= 1;
                else if (next_job[-4] == '\r' && next_job[-3] == '\n')
                        next_job -= 2;
        }
	return next_job;
}
/*get the correct starting point*/
char *get_right_sp (char *start_point ,char type)
{
	start_point = strchr(start_point,'\n')+1;
	return start_point;	
}
/*count useful characters for fasta reads*/
char *fa_count (char *start, int length)
{
	char *reads = (char *) calloc (length+1, sizeof(char));
	char *p = reads;
	// conservatively allocate memory
	while (length>0)
	{
		if (*start!='\n')
		{
			p[0]=start[0];
			p++;
		}
		start++;
		length--;
	} 
	p[0] = '\0';
	return reads;
}
