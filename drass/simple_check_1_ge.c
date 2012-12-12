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
  if (help == 1) {
      check_help ();
      exit (1);
  }
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
      if (tole_rate==0)
          tole_rate = mco_suggestion(bl_2->k_mer);
#pragma omp parallel
      {
#pragma omp single nowait
	{
	  while (head != tail) {
#pragma omp task firstprivate(head)
	      {
		if (head->location!=NULL)
                  {
                  //printf ("location->%0.6s\n",head->location);
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
  float contamination_rate =
    (float) (File_head->reads_contam) / (float) (File_head->reads_num);

// JSON output format by default
  printf("{\n");
  printf ("\t\"total_read_count\": %d,\n", File_head->reads_num);
  printf ("\t\"contaminated_reads\": %d,\n", File_head->reads_contam);
  printf ("\t\"contamination_rate\": %f,\n", contamination_rate);
  printf ("\t\"bloom_filename\":\"%s\"\n", filename);
  printf("}\n");

#ifdef DEBUG
  strcat (detail, "Bloomfile\tAll\tContam\tcontam_rate\n");
  strcat (detail, filename);
#endif

  sprintf (buffer, "  %d\t%d\t%f\n", File_head->reads_num,
	   File_head->reads_contam, contamination_rate);
  strcat (detail, buffer);
}

/*-------------------------------------*/
void
statistic_save (char *detail, char *filename, char *prefix)
{
  char *save_file = NULL;
  save_file = prefix_make (filename, NULL, prefix);
  strcat (save_file,".info");

#ifdef DEBUG
  printf ("filename->%s\n", filename);
#endif

  printf ("Info name->%s\n", save_file);
  write_result (save_file, detail);
}
