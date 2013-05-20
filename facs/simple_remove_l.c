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
#include "remove_l.h"
#include "file_dir.h"

#ifndef __clang__
#include<omp.h>
#endif

char *clean_l, *contam_l;
<<<<<<< HEAD
=======

/*
 *
 * THIS IS REDUNDANT CODE, SHOULD BE MERGED WITH REMOVE
 *
 */
static int
remove_l_usage (void)
{
  fprintf (stderr, "\nUsage: ./facs remove [options]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "\t-b reference Bloom filter to query against\n");
  fprintf (stderr, "\t-q FASTA/FASTQ file containing the query\n");
  fprintf (stderr,
	   "\t-l input list containing all Bloom filters, one per line\n");
  fprintf (stderr, "\t-r single Bloom filter file\n");
  fprintf (stderr, "\t-t threshold value\n");
  return 1;
}

>>>>>>> fd1a93c681ec92d4019f14642c1cc4dd08616bd8

int
remove_l (char *source, char *ref, char *list, char *prefix)
{
/*-------------------------------------*/
  char *position;
  int type = 1;
<<<<<<< HEAD
  float tole_rate = 0;
=======
>>>>>>> fd1a93c681ec92d4019f14642c1cc4dd08616bd8
  char *clean_l2 = clean_l;
  char *contam_l2 = contam_l;
/*-------------------------------------*/
/*-------------------------------------*/
  bloom *bl_2 = NEW (bloom);
  Queue *head = NEW (Queue);
  Queue *tail = NEW (Queue);
  head->location = NULL;
  head->next = tail;
  Queue *head2 = head;

  F_set *File_head = NEW (F_set);
  File_head = make_list (ref, list);
  F_set *File_head2 = File_head;

  position = mmaping (source);
  type = get_parainfo (position, head);
  clean_l = (char *) malloc (strlen (position) * sizeof (char));
  contam_l = (char *) malloc (strlen (position) * sizeof (char));
  while (File_head)
    {
      load_bloom (File_head->filename, bl_2);
      tole_rate = mco_suggestion (bl_2->k_mer);
      printf ("File name->%s", File_head->filename);
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
<<<<<<< HEAD
		      fasta_process_m (bl_2, head, tail, tole_rate, File_head, 1);
		    else
		      fastq_process_m (bl_2, head, tail, tole_rate, File_head, 1);
=======
		      fasta_process_ml (File_head, bl_2, head, tail, clean_l,
					contam_l, tole_rate);
		    else
		      fastq_process_ml (File_head, bl_2, head, tail, clean_l,
					contam_l, tole_rate);
>>>>>>> fd1a93c681ec92d4019f14642c1cc4dd08616bd8
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

<<<<<<< HEAD
  all_save (File_head2, head2, tail, source, clean_l, clean_l2, contam_l, contam_l2, position, type, prefix);
=======
  all_save (File_head2, head2, tail, source, clean_l, clean_l2, contam_l, contam_l2,
	    position, type, prefix);
>>>>>>> fd1a93c681ec92d4019f14642c1cc4dd08616bd8
  munmap (position, strlen (position));

#ifdef DEBUG
  printf ("finish processing...\n");
#endif

  return 0;
}
<<<<<<< HEAD
=======

/*-------------------------------------*/
void
fastq_process_ml (F_set * File_head, bloom * bl, Queue * info, Queue * tail,
		  char *clean_l, char *contam_l, float tole_rate)
{
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
      read_num = count_read (p, next, 2);
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
	fastq_read_check (p, strchr (p, '\n') - p, 'n', bl, tole_rate,
			  File_head);

      if (result == 0)
	{
#pragma omp critical
	  {
	    memcpy (clean_l, temp_start, temp_end - temp_start);
	    clean_l += (temp_end - temp_start);
	    if (*temp_end != '\0')
	      {
		clean_l[0] = '\n';
		clean_l++;
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
fasta_process_ml (F_set * File_head, bloom * bl, Queue * info, Queue * tail,
		  char *clean_l, char *contam_l, float tole_rate)
{
  int read_num = 0, result = 0, countup = 0;
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

      result = fasta_read_check (p, temp, 'n', bl, tole_rate, File_head);
      if (result == 0)
	{
#pragma omp critical
	  {
	    memcpy (clean_l, p, temp - p);
	    clean_l += (temp - p);
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
  //printf ("all->%d\ncontam_l->%d\n", read_num, read_contam_l);
}

/*-------------------------------------*/
void
save_result_ml (char *source, char *obj_file, char *data, char *data2,
		int flag, int type, char *prefix)
{
  printf ("saving...\n");
  char *match = (char *) malloc (4 * HUN * sizeof (char)),
    *so_name = (char *) malloc (2 * HUN * sizeof (char)),
    *obj_name = (char *) malloc (2 * HUN * sizeof (char));

  memset (match, 0, 4 * HUN);
  memset (so_name, 0, 2 * HUN);
  memset (obj_name, 0, 2 * HUN);

  char *so = strrchr (source, '/');
  if (so)
    so += 1;
  //((so = strrchr (source, '/'))) && (so += 1, 1) || (so = NULL);

  char *obj = strrchr (obj_file, '/');
  if (obj)
    obj += 1;
  //((obj = strrchr (obj_file, '/'))) && (obj += 1, 1) || (obj = NULL);

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
    strcat (match, "_clean_l");
  else
    strcat (match, "_contam_l");
  if (type == 1)
    {
      strcat (match, ".fasta");
    }
  else
    {
      strcat (match, ".fastq");
    }
  //printf ("match->%s\n", match);

  write_result (match, data2);

  free (match);

  free (so_name);

  free (obj_name);

  memset (data2, 0, strlen (data2));

  data = data2;

}

>>>>>>> fd1a93c681ec92d4019f14642c1cc4dd08616bd8
/*-------------------------------------*/
int
count_read (char *data, char *next, int type)
{
  //printf ("count_read\n");
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
all_save (F_set * File_head2, Queue * head2, Queue * tail, char *source,
	  char *clean_l, char *clean_l2, char *contam_l, char *contam_l2,
	  char *position, int type, char *prefix)
{
  char *pos, *next, *temp_next;
  int countup;
  Queue *head;
  printf ("clean_l->%s\n", clean_l);
<<<<<<< HEAD
  save_result (source, File_head2->filename, 0, prefix, clean_l, clean_l2, contam_l, contam_l2);	// save the clean_l data
=======
  save_result_ml (source, File_head2->filename, clean_l, clean_l2, 0, type, prefix);	// save the clean_l data
>>>>>>> fd1a93c681ec92d4019f14642c1cc4dd08616bd8
  free (clean_l2);

//File_head2 = File_head2->next;
  printf ("1_dollar_%s\n", File_head2->filename);
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
		  memcpy (contam_l, pos, temp_next - pos);
		  contam_l += temp_next - pos;
		}
	      countup++;
	      pos = temp_next;
	    }
	  head = head->next;
	  //save_result (source,File_head2->filename,contam_l,contam_l2,1);
	}
      memset (clean_l2, 0, strlen (position));
      memset (contam_l2, 0, strlen (position));
<<<<<<< HEAD
      save_result (source, File_head2->filename, 0, prefix, clean_l, clean_l2, contam_l, contam_l2);
=======
      save_result_ml (source, File_head2->filename, contam_l, contam_l2, 1, type,
		   prefix);
>>>>>>> fd1a93c681ec92d4019f14642c1cc4dd08616bd8
      File_head2 = File_head2->next;
    }
}
