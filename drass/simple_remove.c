#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*-------------------------------------*/
//for file mapping in Linux and timing
#include<fcntl.h>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/time.h>
#include<sys/mman.h>
#include<sys/types.h>
/*-------------------------------------*/
#include "tool.h"
#include "bloom.h"
#include "remove.h"
#include "file_dir.h"
/*-------------------------------------*/
//openMP library
#include<omp.h>
//#include<mpi.h>
/*-------------------------------------*/
char *clean, *contam;
/*-------------------------------------*/

int remove_main(int argc, char** argv)
{
  if (argc < 2) remove_help();
/*-------defaults for bloom filter building-------*/ 
  int opt;
  float tole_rate = 0;
  char* ref = NULL;
  char* list = NULL;
  char* target_path = NULL;
  char* source = NULL;
  while ((opt = getopt (argc, argv, "t:r:o:q:l:h")) != -1) {
      switch (opt) {
          case 't':
              (optarg) && ((tole_rate = atof(optarg)), 1);
              break;
          case 'o':    
              (optarg) && ((target_path = optarg), 1);
              break;
          case 'q':  
              (optarg) && (source = optarg, 1);  
              break;
          case 'r':  
              (optarg) && (ref = optarg, 1);  
              break;
          case 'l':
              (optarg) && (list = optarg, 1);  
              break;
          case 'h':
              remove_help();
          case '?':
              printf ("Unknown option: -%c\n", (char) optopt);
              remove_help();
      } 
  } 
  return remove_reads(source, ref, list, target_path, tole_rate);
}
int remove_reads(char *source, char *ref, char *list, char *prefix, float tole_rate)
{
  /*-------------------------------------*/
  int type = 1;
  char *position;
  //char *clean;
  //char *contam;
  char *clean2;
  char *contam2;
  /*-------------------------------------*/
  bloom *bl_2 = NEW (bloom);
  Queue *head = NEW (Queue);
  Queue *tail = NEW (Queue);
  head->next = tail;
  Queue *head2 = head;
  F_set *File_head = NEW (F_set);
  File_head = make_list (ref, list);
  /*-------------------------------------*/
  position = mmaping (source);
  type = get_parainfo (position, head);
  clean = (char *) malloc (strlen (position) * sizeof (char));
  contam = (char *) malloc (strlen (position) * sizeof (char));
  clean2 = clean;
  contam2 = contam;
  /*-------------------------------------*/
  while (File_head)
    {
      memset (clean2, 0, strlen (position));
      memset (contam2, 0, strlen (position));
      load_bloom (File_head->filename, bl_2);
      
      if (tole_rate==0)
      	tole_rate = mco_suggestion (bl_2->k_mer);
#pragma omp parallel
      {
#pragma omp single nowait
	{
	  while (head != tail)
	    {
#pragma omp task firstprivate(head)
	      {
		if (head->location!=NULL)
		  if (type == 1)
		    fasta_process_m (bl_2, head, tail, tole_rate, File_head);
		  else
		    fastq_process_m (bl_2, head, tail, tole_rate, File_head);
	    }
          }
	      head = head->next;
	}			// End of single - no implied barrier (nowait)
      }				// End of parallel region - implied barrier
      save_result (source, File_head->filename, type, prefix, clean, clean2,
		   contam, contam2);
      File_head = File_head->next;
      head = head2;
      bloom_destroy (bl_2);
    }				//end while
  munmap (position, strlen (position));
  printf ("finish processing...\n");
  return 0;
}

/*-------------------------------------*/
void
fastq_process_m (bloom * bl, Queue * info, Queue * tail, float tole_rate, F_set *File_head)
{

  int read_num = 0, read_contam = 0;
  char *p = info->location;
  char *next, *temp_start, *temp_end, *temp_piece = NULL;

  if (info->next == NULL)
    return;

  else if (info->next != tail)
    next = info->next->location;

  else
    next = strchr (p, '\0');

  while (p != next)
    {

      read_num++;

      temp_start = p;

      if (p == '\0' || p == NULL)
	break;

      p = strchr (p, '\n') + 1;

      temp_end = strstr (p, "\n@");

      if (!temp_end)
	temp_end = strchr (p, '\0');
      int result =
	fastq_read_check (p, strchr (p, '\n') - p, 'n', bl, tole_rate, File_head);

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
      else if (result > 0)
	{
#pragma omp critical
	  {
            read_contam++;
	    memcpy (contam, temp_start, temp_end - temp_start);
	    contam += (temp_end - temp_start);
	    if (*temp_end != '\0')
	      {
		contam[0] = '\n';
		contam++;
	      }
	  }
	}


      if (*temp_end == '\0')
	break;

      p = temp_end + 1;

    }				// outside while
//free(key);
  if (temp_piece)
    free (temp_piece);
}

/*-------------------------------------*/
void
fasta_process_m (bloom * bl, Queue * info, Queue * tail, float tole_rate, F_set *File_head)
{
  printf ("fasta processing...\n");

  int read_num = 0, read_contam = 0;

  char *p = info->location;

  char *next;

  char *temp = p;

  if (info->next == NULL)
    return;
  else if (info->next != tail)
    next = info->next->location;
  else
    next = strchr (p, '\0');

  while (p != next)
    {
      read_num++;
      temp = strchr (p + 1, '>');
      if (!temp)
	temp = next;

      int result = fasta_read_check (p, temp, 'n', bl, tole_rate, File_head);
      if (result == 0)
	{
#pragma omp critical
	  {
	    memcpy (clean, p, temp - p);
	    clean += (temp - p);
	  }
	}
      else if (result > 0)
	{
#pragma omp critical
	  {
            read_contam++;
	    memcpy (contam, p, temp - p);
	    contam += (temp - p);
	  }
	}
      p = temp;
    }
  printf ("all->%d\ncontam->%d\n", read_num, read_contam);
}

/*-------------------------------------*/
void
save_result (char *source, char *obj_file, int type, char *prefix,
	     char *clean, char *clean2, char *contam, char *contam2)
{
  printf ("source->%s\n",source);
  printf ("obj_file->%s\n",obj_file);
  printf ("prefix->%s\n",prefix);
  char *so = NULL, *obj = NULL;

  char *match = (char *) malloc (400 * sizeof (char)),
    *mismatch = (char *) malloc (400 * sizeof (char)),
    *so_name = (char *) malloc (200 * sizeof (char)),
    *obj_name = (char *) malloc (200 * sizeof (char));

  memset (match, 0, 400);
  memset (mismatch, 0, 400);
  memset (so_name, 0, 200);
  memset (obj_name, 0, 200);
  
  so = strrchr (source, '/');
  obj = strrchr (obj_file,'/');
  if (so)
     so += 1;
  else 
     so = NULL;
  if (obj)
     obj += 1;
  else
     obj = NULL;
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
      strcat (mismatch, prefix);
    }
  else if (so)
    {
      strncat (match, source, so - source);
      strncat (mismatch, source, so - source);
    }
  //printf ("objname->%s\n",obj_name);
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  strcat (match, so_name);
  strcat (mismatch, so_name);
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  strcat (match, "_");
  strcat (mismatch, "_");
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  strcat (match, obj_name);
  strcat (mismatch, obj_name);
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  strcat (match, "_contam");
  strcat (mismatch, "_clean");

  if (type == 1)
    {
      strcat (match, ".fasta");
      strcat (mismatch, ".fasta");
    }
  else
    {
      strcat (match, ".fastq");
      strcat (mismatch, ".fastq");
    }
  printf ("match->%s\n", match);
  printf ("mis->%s\n", mismatch);

  write_result (match, contam2);
  write_result (mismatch, clean2);
  free (match);
  free (mismatch);
  free (so_name);
  free (obj_name);
  memset (contam2, 0, strlen (contam2));
  memset (clean2, 0, strlen (clean2));
  clean = clean2;
  contam = contam2;
}

/*-------------------------------------*/
