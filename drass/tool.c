#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include "tool.h"
#include "bloom.h"
#include "file_dir.h"
/*-------------------------------------*/

int
fastq_read_check (char *begin, int length, char *model, bloom * bl,
		  float tole_rate)
{
  char *p = begin;

  int distance = length;
  int signal = 0, result = 0;
  char *previous, *key = (char *) malloc (bl->k_mer * sizeof (char) + 1);

#ifdef DEBUG
  printf ("fastq read check...\n");
#endif

  while (distance > bl->k_mer) {
      if (signal == 1) break;

      if (distance >= bl->k_mer) {
          memcpy (key, p, sizeof (char) * bl->k_mer);	//need to be tested
          key[bl->k_mer] = '\0';
          previous = p;
          p += bl->k_mer;
          distance -= bl->k_mer;
  	  } else {
          memcpy (key, previous + distance, sizeof (char) * bl->k_mer);
          p += (bl->k_mer - distance);
          signal = 1;
	  }

#ifdef DEBUG
  printf ("%s\n", key);
#endif

      if (model == "reverse")
	      rev_trans (key);

      if (bloom_check (bl, key)) {
          result = fastq_full_check (bl, begin, length, model, tole_rate);
          if (result > 0)
            return result;
          else if (model == "normal")
            break;
	  }
    } //outside while

  if (model == "reverse")
    return 0;
  else
    return fastq_read_check (begin, length, "reverse", bl, tole_rate);
}

/*-------------------------------------*/
int
fastq_full_check (bloom * bl, char *p, int distance, char *model,
		  float tole_rate)
{

#ifdef DEBUG
  //printf ("fastq full check...\n");
#endif

  int length = distance;
  int count = 0, match_s = 0, mark = 1;
  float result;
  char *previous, *key = (char *) malloc (bl->k_mer * sizeof (char) + 1);

  while (distance >= bl->k_mer)
    {
      memcpy (key, p, sizeof (char) * bl->k_mer);
      key[bl->k_mer] = '\0';
      previous = p;
      p += 1;

      if (model == "reverse")
	rev_trans (key);

      if (count >= bl->k_mer)
	{
	  mark = 1;
		  count = 0;
		}
	      if (bloom_check (bl, key))
		{
		  //printf("hit--->\n");
		  if (mark == 1)
		    {
		      match_s += (bl->k_mer - 1);
		      mark = 0;
		    }
		  else
		    match_s++;
		}
	      else
		//printf("unhit-->\n");
		distance--;
	    }				// end while
	  free (key);
	  result = (float) match_s / (float) length;
	  if (result >= tole_rate)
	    return match_s;
	  else
	    return 0;
	}

	/*-------------------------------------*/
	int
	fasta_read_check (char *begin, char *next, char *model, bloom * bl,
			  float tole_rate)
	{

	  char *p = strchr (begin + 1, '\n') + 1;

	  if (!p || *p == '>')
	    return 1;

	  char *start = p;
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

	      if (model == "reverse")
		rev_trans (key);

	      if (bloom_check (bl, key))
		{
		  result = fasta_full_check (bl, begin, next, model, tole_rate);
		  if (result > 0)
		    return result;
		  //else if (model == "normal")     //use recursion to check the sequence forward and backward
		  //    return fasta_read_check (begin, next, "reverse", bl);
		  else if (model == "normal")
		    break;
		}

	      //memset (key, 0, bl->k_mer);
	    }				//outside while
	  if (model == "reverse")
	    return 0;
	  else
	    return fasta_read_check (begin, next, "reverse", bl, tole_rate);
	}

	/*-------------------------------------*/
	int
	fasta_full_check (bloom * bl, char *begin, char *next, char *model,
			  float tole_rate)
	{
	  int match_s = 0, count = 0, mark = 1;

	  int n = 0, m = 0, count_enter = 0;

	  float result = 0;

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

	      if (model == "reverse")
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
		      //printf("hit--->\n");
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
	  result = (float) match_s / (float) (next - begin - count_enter);
	  if (result >= tole_rate)	//match >tole_rate considered as contaminated
	    return match_s;
	  else
	    return 0;
	}

	/*-------------------------------------*/
	int
	get_parainfo (char *full, Queue * head)
	{
	  printf ("distributing...\n");
	  int type;
          char *previous = NULL;
	  char *temp = full;
	  int cores = omp_get_num_procs ();
	  int offsett = strlen (full) / cores;
	  int add = 0;

	  //printf ("task->%d\n", offsett);

	  Queue *pos = head;

	  if (*full == '>')
	    type = 1;
	  else if (*full == '@')
	    type = 2;
	  else
	    {
	      perror ("wrong format\n");
	      exit (-1);
	    }
	  if (type == 1)
	    {
	      for (add = 0; add < cores; add++)
		{
		  Queue *x = NEW (Queue);

		  if (add == 0 && *full != '>')

		    temp = strchr (full, '>');	//drop the possible fragment

		  if (add != 0)
		    temp = strchr (full + offsett * add, '>');
                  //if (temp)
                  //if (previous!=temp)
                  //    previous = temp;

                  //if (previous!=temp)
                  //{
		  x->location = temp;
		  x->number = add;
		  x->next = pos->next;
		  pos->next = x;
		  pos = pos->next;
                  //}
		}
	    }
	  else
	    {

	      for (add = 0; add < cores; add++)
		{
		  Queue *x = NEW (Queue);
		  x->location = NULL;
		  if (add == 0 && *full != '@')
		    temp = strstr (full, "\n@") + 1;	//drop the fragment

		  //printf("offset->%d\n",offsett*add);
		  
		  if (add != 0)
                    {
		    temp = fastq_relocate(full,offsett*add);
                    /*
                    if (temp)
                        {
			temp++;
                        //if (previous!=temp)
                        //previous = temp;
                        }
                    */
                    } 
                   
          if (previous!=temp)
          {
          previous = temp;
	  x->location = temp;
          //printf ("task->%0.6s\n",x->location);
	  x->number = add;
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
	point = strstr (target + 1, "\n@") + 1;	//point to @

      if (point)
	target = point;
    }
  return target;
}

/*-------------------------------------*/
char *fastq_relocate (char *data, int offset){
     char *target=NULL;
     target = strstr (data + offset, "\n+");
     if (!target)
         return NULL;
     else
         {
         target = strchr (target+1,'\n')+1; 
         if (target!=NULL)
             target = strchr (target+1,'\n')+1;
         
         }
     
     return target;
}
