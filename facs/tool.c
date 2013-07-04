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
/*quick pass for fastq reads using k-mer and 0 overlap*/
int fastq_read_check (char *begin, int length, char mode, bloom * bl, float tole_rate, F_set * File_head)
{
	//printf("mode->%c---read_check\n",mode);
	if (mode == 'r') // make a copy of the read for reverse compliment process
	{
		char *re_compliment = (char *) malloc (sizeof (char) *(length+1));
		re_compliment[length]='\0';
		memcpy(re_compliment, begin, length);
		begin = re_compliment;
		rev_trans (begin,length);
		//printf("reverse->%s\n",begin);
	}
	// initialization
	int result = 0, read_length = length;
	char *start_point = begin;
	while (read_length > 0)
	{
		if (read_length >= bl->k_mer)
		{
			//start_point += bl->k_mer;
			read_length -= bl->k_mer;
		}
      		else
		{
	  		start_point -= (bl->k_mer-read_length);
			read_length = 0;
		}
		if (bloom_check (bl, start_point))
		{
			result = fastq_full_check (bl, begin, length, tole_rate, File_head);
	  		if (result > 0)
			{	
				if (mode == 'r') //free the reverse compliment read copy
					free(begin);
	    			return result;
			}
	  		else if (mode == 'n')
	    			break;
		}
		start_point+=bl->k_mer;
	}		//outside while
  	if (mode == 'r')
	{
		free(begin); //free the reverse compliment read copy
    		return 0;
	}	
  	else
    		return fastq_read_check (begin, length, 'r', bl, tole_rate, File_head);
}
/*full check for fastq sequence with k-mer and k-1 overlap*/
int fastq_full_check (bloom * bl, char *start_point, int length, float tole_rate, F_set * File_head)
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
			//printf("%0.15s\n",start_point);
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
	result = (float) (match_time * bl->k_mer + conse) / (float) (length * bl->k_mer - 2 * bl->dx + length - bl->k_mer + 1);
	#pragma omp atomic
	File_head->hits += match_time;
	#pragma omp atomic
	File_head->all_k += (length - bl->k_mer);
	if (result >= tole_rate)
		return match_s;
	else
		return 0;
}
/*fasta read quick check using k-mer and 0 overlap*/
int
fasta_read_check (char *begin, char *next, char model, bloom * bl, float tole_rate, F_set * File_head)
{

  char *p = strchr (begin + 1, '\n') + 1;

  if (!p || *p == '>')
    return 1;

  int n, m, result, count_enter=0;
  char *key = (char *) malloc ((bl->k_mer + 1) * sizeof (char));
  char *pre_key = (char *) malloc ((bl->k_mer + 1) * sizeof (char));

  key[bl->k_mer] = '\0';

  while (p != next)
    {
      while (n < bl->k_mer)
	{
	  if (p[m] == '>' || p[m] == '\0')
	    {
	      m--;
	      break;
	    }

	  if (p[m] != '\r' && p[m] != '\n')
	    key[n++] = p[m];
	  else
	    count_enter++;
	  m++;
	}			//inner while

      if (m == 0)
	break;

      if (strlen (key) == bl->k_mer)
	memcpy (pre_key, key, sizeof (char) * (bl->k_mer + 1));

      else
	{
	  char *temp_key = (char *) malloc (bl->k_mer * sizeof (char));

	  memcpy (temp_key, pre_key + strlen (key), bl->k_mer - strlen (key));

	  memcpy (temp_key + bl->k_mer - strlen (key), key, sizeof (char) * (strlen (key) + 1));

	  free (key);

	  key = temp_key;

	}
      p += m;

      n = 0;

      m = 0;

      //if (model == 'r')
	//rev_trans (key);

      if (bloom_check (bl, key))
	{
	  result = fasta_full_check (bl, begin, next, model, tole_rate, File_head);
	  if (result > 0)
	    return result;
	  //else if (model == 'n')     //use recursion to check the sequence forward and backward
	  //    return fasta_read_check (begin, next, 'r', bl);
	  else if (model == 'n')
	    break;
	}

      //memset (key, 0, bl->k_mer);
    }				//outside while
  if (model == 'r')
    return 0;
  else
    return fasta_read_check (begin, next, 'r', bl, tole_rate, File_head);
}

/*fasta full check using k-mer and k-1 overlap*/
int
fasta_full_check (bloom * bl, char *begin, char *next, char model, float tole_rate, F_set * File_head)
{
  int match_s = 0, count = 0, mark = 1;

  int n = 0, m = 0, count_enter = 0, match_time = 0;

  short previous = 0, conse = 0;

  float result;

  char *key = (char *) malloc ((bl->k_mer + 1) * sizeof (char));

  begin = strchr (begin + 1, '\n') + 1;

  char *p = begin;

  while (p != next)
    {
      if (*p == '\n')
	count_enter++;
      p++;
    }

  p = begin;

  while (*p != '>' && *p != '\0')
    {
      while (n < bl->k_mer)
	{
	  if (p[m] == '>' || p[m] == '\0')
	    {
	      m--;
	      break;
	    }

	  if (p[m] != '\r' && p[m] != '\n')
	    key[n++] = p[m];

	  m++;
	}
      key[n] = '\0';

     // if (model == 'r')
	//rev_trans (key);
      //printf("key->%s\n",key);
      if (count >= bl->k_mer)
	{
	  mark = 1;
	  count = 0;
	}
      //if (strlen (key) == bl->k_mer)
	//{
	  if (bloom_check (bl, key))
	    {
	      match_time++;
	      if (previous == 1)
		conse++;
	      else
		{
		  conse += bl->k_mer;
		  previous = 1;
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
	      previous = 0;
	      //printf("unhit--->\n");
	    }

	  count++;
	//}			//outside if
      //printf("score->%d\n",match_s);
      p++;
      if (p[0] == '\n')
	p++;
      n = 0;
      m = 0;
    }				// end of while
  //result = (float) match_s / (float) (next - begin - count_enter);
  //result = (float) match_time*(bl->k_mer)/(float)((next-begin-count_enter-bl->k_mer+2)*(bl->k_mer)+2*dx_add(bl->k_mer));
  //result = (float) ((match_time+conse)*(bl->k_mer))/(float)((next-begin-count_enter-bl->k_mer+2+conse)*(bl->k_mer)+2*dx_add(bl->k_mer));
  //result = (float) ((match_time)*(bl->k_mer))/(float)((next-begin-count_enter-bl->k_mer+2)*(bl->k_mer)+2*dx_add(bl->k_mer));
  //result = (float)(match_time*bl->k_mer+conse)/(float)((next-begin-count_enter-bl->k_mer+2)*bl->k_mer+conse+2*dx_add(bl->k_mer));
  //printf ("result1->%f\n",result);
  //result = (float)(match_time*bl->k_mer)/(float)((next-begin-count_enter)*bl->k_mer-2*dx_add(bl->k_mer-1));
  result = (float) (match_time * bl->k_mer + conse) / (float) ((next - begin - count_enter) * bl->k_mer - 2 * bl->dx + (next - begin - count_enter) - bl->k_mer + 1);

#pragma omp atomic
  File_head->hits += match_time;
#pragma omp atomic
  File_head->all_k += (next - begin - count_enter - bl->k_mer);

  if (result >= tole_rate)	//match >tole_rate considered as contaminated
    return match_s;
  else
    return 0;
}

/*Parallel job distribution*/
int
get_parainfo (char *full, Queue * head)
{
#ifdef DEBUG
  printf ("distributing...\n");
#endif
	  int type = 0;
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
       //   Queue *x = NEW (Queue);
      cores = 1;
      int length = 0;
      if (full != NULL) {
          offset = strlen(full) / cores;
          if (*full == '>')
            type = 1;
          else if (*full == '@')
            type = 2;
          else {
                fprintf(stderr, "File format not supported\n");
                exit(EXIT_FAILURE);
          }
      }
      
      if (type == 1) {
              for (add = 0; add < cores; add++) {
                  Queue *x = NEW (Queue);
                  if (add == 0 && *full != '>')
                    temp = strchr (full, '>');	//drop the possible fragment

                  if (add != 0)
                    temp = strchr (full + offset * add, '>');
                  x->location = temp;
                  x->number = &add;
                  x->next = pos->next;
                  pos->next = x;
                  pos = pos->next;
              }

	  } else {
              //char *tx = strchr(full,'\n');
              //length = strchr(tx+1,'\n')-(tx+1);
     
	      for (add = 0; add < cores; add++) {
              char *tx = strchr(full,'\n');
              length = strchr(tx+1,'\n')-(tx+1);
              
	      Queue *x = NEW (Queue);
              x->location = NULL;
              //char *tx = strchr(full,'\n');
              //length = strchr(tx+1,'\n')-(tx+1);
              if (add != 0)
                  temp = fastq_relocate(full, offset*add, length);
                       
              if (previous!=temp) {
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
char *
jump (char *target, int type, float sampling_rate)
{
  //printf("here\n");
  float seed = rand () % 10;

  if (seed >= (float) sampling_rate * 10)
    {

      char *point;

      if (type == 1)
	point = strchr (target + 1, '>');	//point to >
      else
	{
	  point = strstr (target + 1, "\n+") + 1;	//point to +
	  point = strchr (point, '\n') + 1;	//point to quality line
	  point = strchr (point, '\n') + 1;	//point to next read starting
	}
      if (point)
	target = point;
    }
  return target;
}


/*relocate the starting points (correct @ positions) for fastq files*/
char *
fastq_relocate (char *data, int offset, int length)
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


/*scoring system scheme*/
int
dx_add (int k_mer)
{
  int x;
  int y = 0;
  for (x = 1; x < k_mer; x++)
    y += x;
  return y;
}
/*get read length for fastq file*/
int
fq_read_length (char *data)
{
  char *origin = data;
  while (*data != '\n')
    data--;
  return origin - data;
}
