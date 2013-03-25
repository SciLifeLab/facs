#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libgen.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>

#include "tool.h"
#include "check.h"
#include "bloom.h"
#include "file_dir.h"

#include <omp.h>

void
fastq_process (bloom * bl, Queue * info, Queue * tail, F_set * File_head,
	       float sampling_rate, float tole_rate)
{

  char *p = info->location;
  char *next = NULL, *temp = NULL, *temp_piece = NULL;

  if (info->location[0] != '@')
    {
      return;
    }
  else if (info->next != tail && info->next->location != NULL)
    {
      next = info->next->location;
    }
  else
    {
      next = strchr (p, '\0');
      //if (next[-2]=='\r')
      //printf ("%d\n",next[-1]);
      if ((next[-1]) == '\n' && next[-2] == '\n')
	next -= 1;
      else if (next[-4] == '\r' && next[-3] == '\n')
	next -= 2;
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
      if (fastq_read_check
	  (p, strchr (p, '\n') - p, 'n', bl, tole_rate, File_head) > 0)
	{
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
  char *temp_next, *next, *temp;

  if (info->location == NULL)
    return;
  else if (info->next != tail)
    next = info->next->location;
  else
    {
      next = strchr (info->location, '\0');
      //if ((next-1)=='\n')
      //next -= 1;
      if ((next[-1]) == '\n' && next[-2] == '\n')
	next -= 1;
      else if (next[-4] == '\r' && next[-3] == '\n')
	next -= 2;
    }

  char *p = info->location;

  while (p != next)
    {
      temp = jump (p, 1, sampling_rate);	// generate random number and judge
                                            // if need to scan this read

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

      if (fasta_read_check (p, temp_next, 'n', bl, tole_rate, File_head) > 0)
	{
#pragma omp atomic
	  File_head->reads_contam++;
	}

      p = temp_next;
    }
}

void
evaluate (char *detail, char *filename, F_set * File_head, char* query)
{
  char buffer[200] = { 0 };
  float contamination_rate = (float) (File_head->reads_contam) /
                             (float) (File_head->reads_num);

// JSON output format by default
  printf ("{\n");
  printf ("\"%s\": ", basename(query)); //sample (query)
  printf ("{\n");
  printf ("\"%s\": ", basename(filename)); //reference
  printf ("{\n");
  printf ("\t\"total_read_count\": %lld,\n", File_head->reads_num);
  printf ("\t\"contaminated_reads\": %lld,\n", File_head->reads_contam);
  printf ("\t\"total_hits\": %lld,\n", File_head->hits);
  printf ("\t\"contamination_rate\": %f,\n", contamination_rate);
  printf ("\t\"bloom_file\":\"%s\"\n", filename);
  printf ("}\n");
  printf ("}\n");
  printf ("}\n");

#ifdef DEBUG
  strcat (detail, "Bloomfile\tAll\tContam\tcontam_rate\n");
  strcat (detail, filename);
#endif

  sprintf(buffer, "  %lld\t%lld\t%f\n", File_head->reads_num,
          File_head->reads_contam, contamination_rate);
  strcat (detail, buffer);
}

void
statistic_save (char *detail, char *filename, char *prefix)
{
  char *save_file = NULL;
  int length = 0;
  //printf ("prefix->%s\n",prefix);
  if (prefix!=NULL && prefix[0]=='.') {
      prefix+=2;
      length = strrchr(prefix,'/')-prefix+1;
      if (length != 0 && strrchr(prefix,'/')!=NULL) {
           save_file =(char *) calloc (length, sizeof (char));
           memcpy(save_file,prefix,length);
           prefix = save_file;
           save_file = NULL;
	  } else {
           prefix = NULL;
	  }
  }
  if (prefix!=NULL)
      if (prefix[strlen(prefix)-1]=='/')
          prefix[strlen(prefix)-1]='\0'; 

  save_file = prefix_make (filename, NULL, prefix);
  
  if (is_dir(prefix) || prefix==NULL)
      strcat (save_file, ".info");
  if (strrchr(save_file,'/')==save_file)
      save_file++;

#ifdef DEBUG
  printf ("Basename->%s\n", filename);
  printf ("Info name->%s\n", save_file);
#endif
  write_result (save_file, detail);
}
