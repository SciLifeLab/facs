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
#include "check.h"
#include "bloom.h"
#include "file_dir.h"
/*-------------------------------------*/
//openMP library
#include<omp.h>
//#include<mpi.h>
/*-------------------------------------*/
int check_main (char *source, char *ref, float tole_rate, float sampling_rate, char *list, char *prefix, int help)
{
  if (help == 1)
    {
      check_help ();
      exit (1);
    }
  /*-------------------------------------*/
  long sec, usec, i;
  
  /*-------------------------------------*/
  char *position;
  char *detail = (char *) malloc (1000 * 1000 * sizeof (char));
  memset (detail, 0, 1000 * 1000);
  int type = 0;
  /*-------------------------------------*/
  Queue *head = NEW (Queue);
  Queue *tail = NEW (Queue);
  bloom *bl_2 = NEW (bloom);
  Queue *head2;
  head->location=NULL;
  head2 = head;
  head->next = tail;

  F_set *File_head = NEW (F_set);
  File_head = make_list (ref, list);
  /*-------------------------------------*/
  position = mmaping (source);
 
  type = get_parainfo (position,head);
  /*-------------------------------------*/
  while (File_head)
    {
      load_bloom (File_head->filename, bl_2);
#pragma omp parallel
      {
#pragma omp single nowait
	{
	  while (head != tail) {
//#pragma omp task firstprivate(head)
	      {
		if (head->location!=NULL)
                  {
		  if (type == 1)
		    fasta_process (bl_2, head, tail, File_head, sampling_rate,
				   tole_rate);
		  else
		    fastq_process (bl_2, head, tail, File_head, sampling_rate,
		  		   tole_rate);
		  }
	      }
	      head = head->next;
	    }
	}			// End of single - no implied barrier (nowait)
      }				// End of parallel region - implied barrier
      evaluate (detail, File_head->filename, File_head);
      /*-------------------------------------*/
      File_head = File_head->next;
      head = head2;
      bloom_destroy (bl_2);
    }				//end while

  statistic_save (detail, source, prefix);
  munmap (position, strlen (position));

  //check ("test.fna","k_12.bloom","r", prefix, 1, 0.8);

  return 0;
}

/*-------------------------------------*/
void
fastq_process (bloom * bl, Queue * info, Queue *tail, F_set * File_head,
	       float sampling_rate, float tole_rate)
{

#ifdef DEBUG
  printf ("fastq processing...\n");
#endif
  
  char *p = info->location;
  char *next, *temp, *temp_piece = NULL;

  if (info->location[0] != '@') {
    return;
  } else if (info->next != tail && info->next->location!=NULL) {
    next = info->next->location;
  } else {
    next = strchr (p, '\0');
  }

  while (p != next)
    {
      temp = jump (p, 2, sampling_rate);	//generate random number and judge if need to scan this read

      if (p != temp)
	{
	  p = temp;
	  continue;
	}

#pragma omp atomic
      File_head->reads_num++;

      p = strchr (p, '\n') + 1;
      if (fastq_read_check (p, strchr (p, '\n') - p, "normal", bl, tole_rate)> 0) {
#pragma omp atomic
	File_head->reads_contam++;
      }

      p = strchr (p, '\n') + 1;
      p = strchr (p, '\n') + 1;
      p = strchr (p, '\n') + 1;
    }				// outside while
  if (temp_piece)
    free (temp_piece);

}

/*-------------------------------------*/
void
fasta_process (bloom * bl, Queue * info, Queue * tail, F_set * File_head,
	       float sampling_rate, float tole_rate)
{

#ifdef DEBUG
  printf ("fasta processing...\n");
#endif
  char *temp_next, *next, *temp, *temp_piece = NULL;

  if (info->location == NULL)
    return;
  else if (info->next != tail)
    next = info->next->location;
  else
    next = strchr (info->location, '\0');

  char *p = info->location;

  while (p != next)
    {
      temp = jump (p, 1, sampling_rate);	//generate random number and judge if need to scan this read

      if (p != temp)
	{
	  p = temp;
	  continue;
	}

#pragma omp atomic
      File_head->reads_num++;

      temp_next = strchr (p + 1, '>');
      if (!temp_next)
	temp_next = next;

      if (fasta_read_check (p, temp_next, "normal", bl, tole_rate) > 0)
	{
#pragma omp atomic
	  File_head->reads_contam++;
	}

      p = temp_next;
    }
}

/*-------------------------------------*/
void
evaluate (char *detail, char *filename, F_set * File_head)
{
  char buffer[200] = { 0 };

  printf ("all->%d\n", File_head->reads_num);
  printf ("contam->%d\n", File_head->reads_contam);
  printf ("bloomname->%s\n", filename);

  float contamination_rate =
    (float) (File_head->reads_contam) / (float) (File_head->reads_num);

  if (contamination_rate == 0)
    printf ("clean data...\n");
  else
    printf ("contamination rate->%f\n", contamination_rate);

  strcat (detail, "\nxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  strcat (detail, "Bloomfile\tAll\tContam\tcontam_rate\n");

  strcat (detail, filename);

  sprintf (buffer, "  %d\t%d\t%f\n", File_head->reads_num,
	   File_head->reads_contam, contamination_rate);
  strcat (detail, buffer);
}

/*-------------------------------------*/
void
statistic_save (char *detail, char *filename, char *prefix)
{
  char *position1,
    *save_file = (char *) malloc (200 * sizeof (char)),
    *possible_prefix = (char *) malloc (100 * sizeof (char));

  memset (save_file, 0, 200);
  memset (possible_prefix, 0, 100);

  position1 = strrchr (filename, '/');

  printf ("filename->%s\n", filename);

  if (!prefix)
    {
      if (position1)
	strncat (possible_prefix, filename, position1 + 1 - filename);
    }
  else
    strcat (possible_prefix, prefix);

  strcat (save_file, possible_prefix);

  if (position1)
    strncat (save_file, position1 + 1, strrchr (filename, '.') - position1);
  else
    strncat (save_file, filename, strrchr (filename, '.') - filename + 1);

  strcat (save_file, "info");

  printf ("Info name->%s\n", save_file);

  write_result (save_file, detail);
}

/*-------------------------------------*/
/*
int
checky (char *query, char *reference, char l, char *target_path,
	double sampling_rate, double tole_rate)
{
  int type;
  char *position;
  Queue *head = NEW (Queue);
  bloom *bl_2 = NEW (bloom);
  Queue *tail, *head2;
  head2 = head;
  head->next = tail;
  F_set *File_head = NEW (F_set);

  if (l == 'l')
    File_head = make_list (NULL, reference);
  else
    File_head = make_list (reference, NULL);
  printf ("File_head->%s\n", File_head->filename);
  char *detail = (char *) malloc (1000 * 1000 * sizeof (char));

  if (strstr (query, ".fifo"))
    position = large_load (query);
  else
    position = mmaping (query);

  type = get_parainfo (position, head);

  while (File_head)
    {
      File_head->reads_num = 0;
      File_head->reads_contam = 0;
      load_bloom (File_head->filename, bl_2);

#pragma omp parallel
      {
#pragma omp single nowait
	{
	  while (head != tail)
	    {
#pragma omp task firstprivate(head)
	      {
		if (head->location)
		  if (type == 1)
		    fasta_process (bl_2, head, tail, File_head, sampling_rate,
				   tole_rate);
		  else
		    fastq_process (bl_2, head, tail, File_head, sampling_rate,
				   tole_rate);

	      }
	      head = head->next;
	    }
	}			// End of single - no implied barrier (nowait)
      }				// End of parallel region - implied barrier
      printf ("num->%d\ncontam->%d\n", File_head->reads_num,
	      File_head->reads_contam);
      evaluate (detail, File_head->filename, File_head);

      File_head = File_head->next;

      head = head2;

      bloom_destroy (bl_2);
    }				//end while

  statistic_save (detail, query, target_path);

  if (!strstr (query, ".fifo"))
    munmap (position, strlen (position));

  return 0;
}
*/
/*-------------------------------------*/
/*
int
check (char *query, char *reference, char l, char *target_path,
       double sampling_rate, double tole_rate)
{
  checky (query, reference, l, target_path, sampling_rate, tole_rate);
  return 0;
}
*/
