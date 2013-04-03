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

#ifndef __clang__ 
#include <omp.h>
#endif

void
fastq_process (bloom * bl, Queue * info, Queue * tail, F_set * File_head,
	       float sampling_rate, float tole_rate)
{

  char *p = info->location;
  char *next = NULL, *temp = NULL, *temp_piece = NULL;

  if (info->location[0] != '@') {
      return;
  } else if (info->next != tail && info->next->location != NULL) {
      next = info->next->location;
  } else {
      next = strchr (p, '\0');
      //if (next[-2]=='\r')
      //printf ("%d\n",next[-1]);
      if ((next[-1]) == '\n' && next[-2] == '\n')
	next -= 1;
      else if (next[-4] == '\r' && next[-3] == '\n')
	next -= 2;
  }

  while (p != next) {
      temp = jump (p, 2, sampling_rate);	//generate random number and judge if need to scan this read

      if (p != temp) {
          p = temp;
          continue;
	  }
#pragma omp atomic
      File_head->reads_num++;

      p = strchr (p, '\n') + 1;
      if (fastq_read_check(p, strchr(p, '\n')-p, 'n',
          bl, tole_rate, File_head) > 0) {
#pragma omp atomic
	      File_head->reads_contam++;
	  }

      p = strchr (p, '\n') + 1;
      p = strchr (p, '\n') + 1;
      p = strchr (p, '\n') + 1;
    }				// outside while

  if (temp_piece)
    free(temp_piece);
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
  else {
      next = strchr (info->location, '\0');
      //if ((next-1)=='\n')
      //next -= 1;
      if ((next[-1]) == '\n' && next[-2] == '\n')
	next -= 1;
      else if (next[-4] == '\r' && next[-3] == '\n')
	next -= 2;
  }

  char *p = info->location;

  while (p != next) {
      temp = jump (p, 1, sampling_rate);	// generate random number and judge
                                            // if need to scan this read
      if (p != temp) {
	    p = temp;
	    continue;
      }

#pragma omp atomic
      File_head->reads_num++;
      temp_next = strchr (p + 1, '>');

      if (!temp_next)
        temp_next = next;

      if (fasta_read_check (p, temp_next, 'n', bl, tole_rate, File_head) > 0) {
#pragma omp atomic
	    File_head->reads_contam++;
	  }

      p = temp_next;
    }
}

void
report(char *detail, char *filename, F_set * File_head, char* query,
       char* fmt, char* prefix)
{
  char buffer[200] = { 0 };
  float contamination_rate = (float) (File_head->reads_contam) /
                             (float) (File_head->reads_num);


  if(!fmt){
      return;
  // JSON output format (via stdout)
  } else if(!strcmp(fmt, "json")) {
      isodate(buffer);

      printf("{\n");
      printf("\t\"timestamp\": \"%s\"\n", buffer);
      printf("\t\"sample\": \"%s\"\n", basename(query)); //sample (query)
      printf("\t\"bloom_filter\": \"%s\"\n", basename(filename)); //reference
      printf("\t\"total_read_count\": %lld,\n", File_head->reads_num);
      printf("\t\"contaminated_reads\": %lld,\n", File_head->reads_contam);
      printf("\t\"total_hits\": %lld,\n", File_head->hits);
      printf("\t\"contamination_rate\": %f,\n", contamination_rate);
      printf("}\n");

  // TSV output format (via file in CWD)
  } else if (!strcmp(fmt, "tsv")) {
      strcat (detail, "sample\tbloom_filter\ttotal_read_count\t\
contaminated_reads\tcontamination_rate\n");

      sprintf(buffer, "%s\t%s\t%lld\t%lld\t%f\n", basename(query),
              basename(filename), File_head->reads_num,
              File_head->reads_contam, contamination_rate);
      strcat(detail, buffer);

      write_result(strcat(basename(query), ".tsv"), detail);
  }
}
