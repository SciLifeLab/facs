#include <zlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>
#include "tool.h"
#include "bloom.h"
#include "check.h"
#include "file_dir.h"
#include "big_query.h"
/*-------------------------------------*/
#include <omp.h>
/*-------------------------------------*/
#define ONEG 1000000000
#define ONE 100
/*-------------------------------------*/

int bq_main(char *source,char *ref,float tole_rate,float sampling_rate,char *list,char *prefix,int help)
{
  /*-------------------------------------*/
  /*
  long sec, usec, i;
  struct timezone tz;
  struct timeval tv, tv2;
  gettimeofday (&tv, &tz);
  */
  /*-------------------------------------*/
  gzFile zip;
  int type = 0;
  BIGCAST offset = 0;
    if (help == 1)
    {
      check_help ();
      exit (1);
    }
  if ((zip = gzopen(source,"rb"))<0)
    {
	 perror("source open error...\n");
         exit(-1);
    }
  char *detail = (char*)malloc((ONE*ONE*ONE)*sizeof(char));
  char *position =  (char*)malloc((ONEG+1)*sizeof(char));
  /*-------------------------------------*/
  if (strstr(source,".fastq")||strstr(source,".fq"))
      type = 2;
  else
      type = 1;
  /*-------------------------------------*/
  Queue *head = NEW (Queue);
  Queue *tail = NEW (Queue);
  bloom *bl_2 = NEW (bloom);
  head->next = tail;
  F_set *File_head = make_list (ref, list);
  File_head->reads_num = 0;
  File_head->reads_contam = 0;
  load_bloom (File_head->filename, bl_2);  //load a bloom filter
  /*-------------------------------------*/
  while (offset!=-1)
    {
	offset = CHUNKer(zip,offset,ONEG,position,type);
	printf("length->%d\n",(int)strlen(position));
  /*-------------------------------------*/
        get_parainfo (position, head);
#pragma omp parallel
{
#pragma omp single nowait
	{
	  while (head != tail)
	    {
#pragma omp task firstprivate(head)
	      {
		if (head->location)
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
  evaluate (detail, File_head->filename, File_head);
      /*-------------------------------------*/
  memset (position, 0, strlen(position));
    }				//end while
  gzclose(zip);
  bloom_destroy (bl_2);
  statistic_save (detail, source, prefix);
  /*
#ifdef DEBUG
  gettimeofday (&tv2, &tz);
  sec = tv2.tv_sec - tv.tv_sec;
  usec = tv2.tv_usec - tv.tv_usec;
#endif
  */
  //printf ("total=%ld sec\n", sec);
  return 0;
}

char *strrstr(char *s, char *str)
{
    char *p; 
    int len = strlen(s);
    for (p = s + len - 1; p >= s; p--) {
        if ((*p == *str) && (memcmp(p, str, strlen(str)) == 0)) 
            return p;
    }   
    return NULL;
}

BIGCAST CHUNKer(gzFile zip,BIGCAST offset,int chunk,char *data,int type)
{
char c, v;
char *pos;
//if (strstr(filename,".fastq")||strstr(filename,".fq"))
if (type == 2)
    v = '@';
else 
    v = '>';

if (offset == 0)
    while (offset <10*ONE)
    {
    c = gzgetc(zip);
    putchar (c);
    	if (c == v)
       	    break;
    offset++;
    }
    
gzseek (zip,offset,SEEK_SET);
gzread (zip,data,chunk);
	  
pos = strrstr (data,"\n@");
offset += (pos-data+1);

if (strlen(data)<chunk)
    offset=-1;
memset(pos+1,0,strlen(pos));

return offset;

}

BIGCAST CHUNKgz(gzFile zip, BIGCAST offset,int chunk,char *position,char *extra,int type)
{
	        memset(position,0,chunk);
  			  char c, *position2 = position;
  			  char *x;
  			  int num=0;
  	      if (offset == 0)
          while (offset <10*ONE)
          {
             c = gzgetc(zip);
    	       if ((c == '@' && type==2)&&(c == '>' && type==1))
       	          break;
             offset++;
          }		
          if (extra!=NULL)
          	  {
          	  memcpy(position,extra,strlen(extra));
          	  position+=strlen(extra);
          	  }
          free(extra);
          while (((c=gzgetc(zip))!=EOF)&&(num<chunk))
          {
             *position=c;
             position++;
             num++;
          }
	        x = strrstr(position2,"\n@");
	        extra = (char*)malloc((position-x+1)*sizeof(char));
	        memcpy (x,extra,position-x+1);
	        offset+=(position-x+1);
return offset;	
}


BIGCAST get_size (char *filename)
{
BIGCAST tim;
struct stat statbuf;
if ((tim=open(filename, O_RDONLY))<0)
    {
     printf("open file error...\n");
     exit(-1);
    }
fstat (tim, &statbuf);
return statbuf.st_size;
}
