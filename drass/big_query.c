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


static int
query_usage(void)
{
    fprintf(stderr, "\nUsage: ./facs query [options]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-r reference bloom filter to query against\n");
    return 1;
}

int bq_main(int argc, char** argv)
{
  if (argc < 2) return query_usage();
  
/*-------defaults for bloom filter building-------*/ 
  int opt;
  float tole_rate = 0;
  float sampling_rate = 1;

  char* prefix = NULL;
  char* ref = NULL;
  char* list = NULL;
  char* target_path = NULL;
  char* source = NULL;

  while ((opt = getopt (argc, argv, "s:t:r:o:q:l:h")) != -1) {
      switch (opt) {
          case 't':
              (optarg) && ((tole_rate = atof(optarg)), 1);
              break;
          case 's':
              (optarg) && ((sampling_rate = atof(optarg)), 1);
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
              return query_usage();
          case '?':
              printf ("Unknown option: -%c\n", (char) optopt);
              return query_usage();
      } 
  } 

  query(source, ref, tole_rate, sampling_rate, list, target_path);
}

int query(char* query, char* bloom_filter, double tole_rate, double sampling_rate, char* list, char* target_path)
{

  gzFile zip;
  int type = 0;
  BIGCAST offset = 0;
  char *detail = (char*) malloc((ONE*ONE*ONE)*sizeof(char));
  char *position =  (char*) malloc((ONEG+1)*sizeof(char));
  
  bloom *bl_2 = NEW (bloom);
  F_set *File_head = make_list (bloom_filter, list);

  File_head->reads_num = 0;
  File_head->reads_contam = 0;
  File_head->filename = bloom_filter;
  load_bloom (File_head->filename, bl_2);  //load a bloom filter
  if (tole_rate==0)
      tole_rate = mco_suggestion (bl_2->k_mer); 
  
  if ((zip = gzopen(query, "rb"))<0) {
	 perror("query open error...\n");
     exit(-1);
  }

  if (strstr(query, ".fastq") || strstr(query, ".fq"))
      type = 2;
  else
      type = 1;
 
  
  while (offset!=-1) {   
	offset = CHUNKer(zip,offset,ONEG,position,type);
    Queue *head = NEW (Queue);
    head->location = NULL;
    Queue *tail = NEW (Queue);
    head->next = tail;
    Queue *head2 = head;
    get_parainfo (position, head);

#pragma omp parallel
{
#pragma omp single nowait
	{
	  while (head != tail)
	    {
#pragma omp task firstprivate(head)
	      {

#ifdef DEBUG
    printf("head->%0.40s\n",head->location);
#endif
		if (head->location) {
            if (type == 1) {
                fasta_process (bl_2, head, tail, File_head, sampling_rate,
                               tole_rate);
		    } else {
                fastq_process (bl_2, head, tail, File_head, sampling_rate,
                               tole_rate);
            }
        }
	  }
	      head = head->next;
	    }       // End of firstprivate
	}			// End of single - no implied barrier (nowait)
}				// End of parallel region - implied barrier
                //evaluate (detail, File_head->filename, File_head);
  
  memset (position, 0, strlen(position));
  clean_list (head2, tail);
    }				//end while

  evaluate (detail, File_head->filename, File_head);
  gzclose(zip);
  bloom_destroy (bl_2);
  statistic_save (detail, query, target_path);
  
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

void clean_list (Queue* head, Queue *tail)
{
Queue *element;
while (head!=tail)
   {
       element = head->next;
       memset(head,0,sizeof(Queue));
       free(head);
       head = element;
   }
free(tail);
}


BIGCAST CHUNKer(gzFile zip,BIGCAST offset,int chunk,char *data,int type)
{
char c, v;
char *pos;
int length = 0;
//if (strstr(filename,".fastq")||strstr(filename,".fq"))
if (type == 2)
    v = '@';
else 
    v = '>';

if (offset == 0)
    while (offset <10*ONE)
    {
	    c = gzgetc(zip);
 //    putchar (c);
    	if (c == v)
       	    break;
    offset++;
    }
    
gzseek (zip,offset,SEEK_SET);
gzread (zip,data,chunk);

if (length>=chunk)
if (type == 2)
    {
	pos = strrstr (data,"\n+");
	pos = bac_2_n (pos-1);
    }
else
    {
	pos = strrchr (data,'>')-1;	
    }

if (pos)
    {
	offset += (pos-data);
        memset (pos,0,strlen(pos));
    }
if (length<chunk)
    offset=-1;

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

char *bac_2_n (char *filename)
{
     while (*filename!='\n')
           filename--;
     filename--;      //move from \n
     while (*filename!='\n')
           filename--;
     filename++;
     return filename;
}
