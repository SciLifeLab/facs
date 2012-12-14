#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*-------------------------------------*/
//for file mapping in Linux
#include<fcntl.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/time.h>
#include<sys/mman.h>
#include<sys/types.h>
/*-------------------------------------*/
#include "tool.h"
#include "bloom.h"
#include "remove_l.h"
#include "file_dir.h"
/*-------------------------------------*/
//openMP library
#include<omp.h>
/*-------------------------------------*/
char *clean, *contam;
/*-------------------------------------*/
int
remove_main_l (float tole_rate, char *source, char *ref, char *list, char *prefix, int help)
{
  if (help == 1)
    {
      remove_l_help ();
      exit (1);
    }
  long sec, usec, i;
  struct timezone tz;
  struct timeval tv, tv2;
  gettimeofday (&tv, &tz);
/*-------------------------------------*/
  char *position;
  int type = 1;
  char *clean2 = clean;
  char *contam2 = contam;
/*-------------------------------------*/
/*-------------------------------------*/
  bloom *bl_2 = NEW (bloom);
  Queue *head = NEW (Queue);
  Queue *tail = NEW (Queue);
  head->next = tail;
  Queue *head2 = head;
  F_set *File_head = NEW (F_set);
  File_head = make_list (ref, list);
  F_set *File_head2 = File_head;

  position = mmaping (source);
  type = get_parainfo (position, head);
  clean = (char *) malloc (strlen (position) * sizeof (char));
  contam = (char *) malloc (strlen (position) * sizeof (char));
  while (File_head)
    {
      load_bloom (File_head->filename, bl_2);
      printf ("File name->%s ", File_head->filename);
      printf ("File number->%d\n", File_head->number);
#pragma omp parallel
      {
#pragma omp single nowait
	{
	  while (head != tail)
	    {
#pragma omp task firstprivate(head)
	      {
		if (head->location)
		  {
		    //printf("location->%0.20s\n",head->location);
		    if (type == 1)
		      fasta_process_ml (File_head, bl_2, head,tail,clean,contam,tole_rate);
		    else
		      fastq_process_ml (File_head, bl_2, head,tail,clean,contam,tole_rate);
		  }
	      }
	      head = head->next;
	    }
	}			// End of single - no implied barrier (nowait)
      }				// End of parallel region - implied barrier
      head = head2;
      bloom_destroy (bl_2);
      File_head = File_head->next;
    }				// End outside while

  all_save (File_head2, head2, tail, source, clean, clean2, contam, contam2,
	    position, type, prefix);
  munmap (position, strlen (position));

  printf ("finish processing...\n");
#ifdef DEBUG
  gettimeofday (&tv2, &tz);
  sec = tv2.tv_sec - tv.tv_sec;
  usec = tv2.tv_usec - tv.tv_usec;
#endif
  printf ("total=%ld sec\n", sec);

  return 0;
}

/*-------------------------------------*/
void
fastq_process_ml (F_set * File_head, bloom * bl, Queue * info,  Queue * tail, char *clean, char *contam, float tole_rate)
{
  printf ("fastq processing...\n");

  int read_num = 0, result = 0, countup = 0;
  char *p = info->location;
  char *next, *temp_start, *temp_end;

  if (info->next == NULL)
    return;
  else if (info->next != tail)
    next = info->next->location;
  else
    next = strchr (p, '\0');

  if (info->score == NULL)
    {
      read_num = count_read (p, next,2);
      printf ("read_num->%d\n", read_num);
      info->score = (short *) malloc (read_num * sizeof (short));
      info->number = (short *) malloc (read_num * sizeof (short));
    }
  while (p != next)
    {
      temp_start = p;

      if (p == '\0' || p == NULL)
	break;

      p = strchr (p, '\n') + 1;

      temp_end = strstr (p, "\n@");

      if (!temp_end)
	temp_end = strchr (p, '\0');

      result =
	fastq_read_check (p, strchr (p, '\n') - p, "normal", bl, tole_rate);

      if (result == 0)
	{
#pragma omp critical
	  {
	    memcpy (clean, temp_start, temp_end - temp_start);
	    clean += (temp_end - temp_start);
	    if (*temp_end != '\0')
	      {
		clean[0] = '\n';
		clean++;
	      }
	  }
	}
      else if (result == -1)
	{
	  continue;
	}

      //(optarg) && (source = optarg, 1);
      if (info->score[countup] < result)
	{
	  info->score[countup] = result;	//record score 
	  info->number[countup] = File_head->number;	//record bloom number
	}
      if (*temp_end == '\0')
	break;
      p = temp_end + 1;
      countup++;
    }				// end while
}

/*-------------------------------------*/
void
fasta_process_ml (F_set * File_head, bloom * bl, Queue * info, Queue * tail, char *clean, char *contam, float tole_rate)
{
  printf ("fasta processing...\n");

  int read_num = 0, result = 0, countup = 0, sign = 0;

  char *p = info->location;

  char *next;

  char *temp = p;

  if (info->next == NULL)
    return;
  else if (info->next != tail)
    next = info->next->location;
  else
    next = strchr (p, '\0');

  if (info->score == NULL)
    {
      read_num = count_read (p, next, 1);
      info->score = (short *) malloc (read_num * sizeof (short));
      info->number = (short *) malloc (read_num * sizeof (short));
    }

  while (p != next)
    {
      temp = strchr (p + 1, '>');
      if (!temp)
	temp = next;

      result = fasta_read_check (p, temp, "normal", bl, tole_rate);
      //printf ("result->%d\n",result);
      if (result == 0)
	{
#pragma omp critical
	  {
	    memcpy (clean, p, temp - p);
	    clean += (temp - p);
	  }
	}
      else if (result == -1)
	continue;
      if (info->score[countup] < result)
	{
	  info->score[countup] = result;	//record score 
	  info->number[countup] = File_head->number;	//record bloom number
	}

      countup++;

      p = temp;
    }				// end while
  //printf ("all->%d\ncontam->%d\n", read_num, read_contam);
}

/*-------------------------------------*/
void
save_result_ml (char *source, char *obj_file, char *data, char *data2,
		int flag, int type, char *prefix)
{
  printf ("saving...\n");
  char *match = (char *) malloc (4*HUN * sizeof (char)),
    *so_name = (char *) malloc (2*HUN * sizeof (char)),
    *obj_name = (char *) malloc (2*HUN * sizeof (char));

  memset (match, 0, 4*HUN);
  memset (so_name, 0, 2*HUN);
  memset (obj_name, 0, 2*HUN);

  char *so;
  ((so = strrchr (source, '/'))) && (so += 1, 1) || (so = NULL);

  char *obj;
  ((obj = strrchr (obj_file, '/'))) && (obj += 1, 1) || (obj = NULL);

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
    }
  else if (so)
    {
      strncat (match, source, so - source);
    }
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  strcat (match, so_name);
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  strcat (match, "_");
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  strcat (match, obj_name);
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  if (flag == 0)
    strcat (match, "_clean");
  else
    strcat (match, "_contam");
  if (type == 1)
    {
      strcat (match, ".fasta");
    }
  else
    {
      strcat (match, ".fastq");
    }
  printf ("match->%s\n", match);

  write_result (match, data2);

  free (match);

  free (so_name);

  free (obj_name);

  memset (data2, 0, strlen (data2));

  data = data2;

}

/*-------------------------------------*/
int
count_read (char *data, char *next, int type)
{
  printf ("count_read\n");
  int number = 1;
  char *pos, *temp_next;
  pos = data;

  while (pos != next)
    {
      if (type == 1)
	temp_next = strchr (pos + 1, '>');
      else
	temp_next = strstr (pos + 1, "\n@");
      number++;
      pos = temp_next;
      if (!pos)
	break;
    }
  return number;
}

/*-------------------------------------*/
void
all_save (F_set * File_head2, Queue * head2, Queue * tail, char *source, char *clean,
	  char *clean2, char *contam, char *contam2, char *position, int type,
	  char *prefix)
{
  char *pos, *next, *temp_next;
  int countup;
  Queue *head;
  save_result_ml (source, File_head2->filename, clean, clean2, 0, type, prefix);	// save the clean data
  free (clean2);

//File_head2 = File_head2->next;
//printf("1_dollar_%s\n",File_head2->filename);
  while (File_head2)
    {
      head = head2;
      head = head->next;
      while (head != tail)
	{
	  countup = 0;
	  pos = head->location;
	  if (head->next->location != NULL)
	    next = head->next->location;
	  else
	    next = strchr (head->location, '\0');
	  //printf("head->%0.20s\n",head->location);                
	  while (pos != next)
	    {
	      if (type == 1)
		temp_next = strchr (pos + 1, '>');
	      else
		temp_next = strstr (pos + 1, "\n@");
	      //printf("temp_next->%0.10s\n",temp_next);
	      //printf("next->%0.10s\n",next);
	      if (temp_next == NULL)
		temp_next = next;
	      //printf("???->%d---%d\n",head->score[countup],head->number[countup]);
	      if (head->score[countup] > 0
		  && (head->number[countup] == File_head2->number))
		{
		  memcpy (contam, pos, temp_next - pos);
		  contam += temp_next - pos;
		}
	      countup++;
	      pos = temp_next;
	    }
	  head = head->next;
	  //save_result (source,File_head2->filename,contam,contam2,1);
	}
      memset (clean2, 0, strlen (position));
      memset (contam2, 0, strlen (position));
      save_result (source, File_head2->filename, contam, contam2, 1, type,
		   prefix);
      File_head2 = File_head2->next;
    }
}
