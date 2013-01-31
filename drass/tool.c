#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <omp.h>
#include "tool.h"
#include "bloom.h"
#include "file_dir.h"
/*-------------------------------------*/

int
fastq_read_check (char *begin, int length, char model, bloom * bl, float tole_rate, F_set *File_head)
{
  char *p = begin;
  int distance = length;
  int signal = 0, result = 0;
  char  *previous, *key = (char *) malloc (bl->k_mer * sizeof (char) + 1);

  while (distance > bl->k_mer)
    {
      if (signal == 1)
	break;

      if (distance >= bl->k_mer)
	{
	  memcpy (key, p, sizeof (char) * bl->k_mer);	//need to be tested
	  key[bl->k_mer] = '\0';
	  p += bl->k_mer;
          previous = p;
	  distance -= bl->k_mer;
	}

      else
	{
	  memcpy (key, previous + distance, sizeof (char) * bl->k_mer);
	  p += (bl->k_mer - distance);
	  signal = 1;
	}

      if (model == 'r')
	rev_trans (key);

      if (bloom_check (bl, key))
	{
	  result = fastq_full_check (bl, begin, length, model, tole_rate, File_head);
	  if (result > 0)
	    return result;
	  else if (model == 'n')
	    break;
	}

    }				//outside while
  if (model == 'r')
    return 0;
  else
    return fastq_read_check (begin, length, 'r', bl, tole_rate, File_head);
}

/*-------------------------------------*/
int
fastq_full_check (bloom * bl, char *p, int distance, char model,
		  float tole_rate, F_set *File_head)
{

  //printf ("fastq full check...\n");

  int length = distance;

  int count = 0, match_s = 0, mark = 1, match_time = 0;

  float result;

  char *key = (char *) malloc (bl->k_mer * sizeof (char) + 1);

  short prev = 0, conse = 0; 

  while (distance >= bl->k_mer)
    {
      memcpy (key, p, sizeof (char) * bl->k_mer);
      key[bl->k_mer] = '\0';
      p += 1;

      if (model == 'r')
	rev_trans (key);

      if (count >= bl->k_mer)
	{
	  mark = 1;
	  count = 0;
	}
 if (strlen (key) == bl->k_mer)
	{
	  if (bloom_check (bl, key))
	    {
	      match_time++;
	      if (prev == 1)
                  conse++;
              else
                  {
                  conse+=bl->k_mer;
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
	      //printf("unhit--->\n");
	    }
	  count++;
	}			//outside if
	distance--;
    }				// end while
  free (key);
  result = (float)(match_time*bl->k_mer+conse)/(float)(length*bl->k_mer-2*bl->dx+conse);
  //result = (float) match_s / (float) length;
  #pragma omp atomic
  File_head->hits+=match_time;
  #pragma omp atomic 
  File_head->all_k+=(length-bl->k_mer);
  if (result >= tole_rate)
    return match_s;
  else
    return 0;
}

/*-------------------------------------*/
int
fasta_read_check (char *begin, char *next, char model, bloom * bl,
		  float tole_rate, F_set *File_head)
{

  char *p = strchr (begin + 1, '\n') + 1;

  if (!p || *p == '>')
    return 1;

  int n, m, result, count_enter;
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

	  memcpy (temp_key + bl->k_mer - strlen (key), key,
		  sizeof (char) * (strlen (key) + 1));

	  free (key);

	  key = temp_key;

	}
      p += m;

      n = 0;

      m = 0;

      if (model == 'r')
	rev_trans (key);

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

/*-------------------------------------*/
int
fasta_full_check (bloom * bl, char *begin, char *next, char model,
		  float tole_rate, F_set *File_head)
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

      if (model == 'r')
	rev_trans (key);
      //printf("key->%s\n",key);
      if (count >= bl->k_mer)
	{
	  mark = 1;
	  count = 0;
	}
      if (strlen (key) == bl->k_mer)
	{
	  if (bloom_check (bl, key))
	    {
	      match_time++;
	      if (previous == 1)
                  conse++;
              else
                  {
                  conse+=bl->k_mer;
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
	}			//outside if
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
  result = (float)(match_time*bl->k_mer+conse)/(float)((next-begin-count_enter)*bl->k_mer-2*bl->dx+conse);
  
  #pragma omp atomic
  File_head->hits+=match_time;
  #pragma omp atomic
  File_head->all_k+=(next-begin-count_enter-bl->k_mer);

  if (result >= tole_rate)	//match >tole_rate considered as contaminated
    return match_s;
  else 
    return 0;
}

int
get_parainfo (char *full, Queue * head)
{
#ifdef DEBUG
	  printf ("distributing...\n");
#endif
	  int type = 0;
          char *previous = NULL;
	  char *temp = full;
	  int cores = omp_get_num_procs ();
	  short add = 0;
          int offset = 0;
	  Queue *pos = head;
       //   Queue *x = NEW (Queue);
          int length = 0;

      if (full != NULL) {
          offset = strlen(full) / cores;
          if (*full == '>')
            type = 1;
          else if (*full == '@')
            type = 2;
          else
            {
            perror ("wrong format\n");
            exit (-1);
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
              char *tx = strchr(full,'\n');
              length = strchr(tx+1,'\n')-(tx+1);
              printf ("reads length->%d\n",length);
	      for (add = 0; add < cores; add++) {
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

/*-------------------------------------*/
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
        point = strchr (point,'\n')+1;          //point to quality line
        point = strchr (point,'\n')+1;          //point to next read starting
        }
      if (point)
	target = point;
    }
  return target;
}

/*-------------------------------------*/
char *fastq_relocate (char *data, int offset, int length){
     char *target=NULL;

     if(data != NULL && offset != 0)
        target = strstr (data + offset, "\n+");

     if (!target)
         return NULL;
     else {
         //if ((strchr(target+1,'\n')-target+1)!=length)
            target = strchr (target+1,'\n')+1; 
         //if (target!=NULL)
            target = strchr (target+1,'\n')+1;
     }
     
     return target;
}
/*-------------------------------------*/
int 
dx_add (int k_mer)
{
   int x;
   int y = 0;
   for (x=1;x<k_mer;x++)
        y+=x;
   return y;
}
