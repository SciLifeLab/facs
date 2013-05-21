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

char *clean_l, *contam_l, *clean_l2, *contam_l2;

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

int
remove_l (char *source, char *ref, char *list, char *prefix)
{
/*-------------------------------------*/
  char *position;
  int type = 1;
  float tole_rate = 0;
  clean_l2 = clean_l;
  contam_l2 = contam_l;
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
		      fasta_process_m (bl_2, head, tail, tole_rate, File_head, 1);
		    else
		      fastq_process_m (bl_2, head, tail, tole_rate, File_head, 1);
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

  all_save (File_head2, head2, tail, source, clean_l, clean_l2, contam_l, contam_l2, position, type, prefix);
  munmap (position, strlen (position));

#ifdef DEBUG
  printf ("finish processing...\n");
#endif

  return 0;
}

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
  save_result (source, File_head2->filename, 0, prefix, clean_l, clean_l2, contam_l, contam_l2);	// save the clean_l data
  free (clean_l2);

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
      save_result (source, File_head2->filename, 0, prefix, clean_l, clean_l2, contam_l, contam_l2);
      File_head2 = File_head2->next;
    }
}
