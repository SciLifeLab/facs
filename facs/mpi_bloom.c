#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64

#include <math.h>
#include <time.h>
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
#include "prob.h"
#include "bloom.h"
#include "big_query.h"
#include "check.h"
#include "hashes.h"
#include "mpi_bloom.h"

/*parallel libraries*/
#include<omp.h>
#include<mpi.h>

/*since the read scanning process is fast enough, a parallel file reading
approach needs to be applied otherwise file reading will become bottleneck.
Back to file mapping approach with multiple access to query file. Zlib has
atomic process with multiple threads reading from a file. So the speed of 
multiple threads reading will be exactly the same as single thread reading
*/
/*-------------------------------------*/
static int mpicheck_usage (void)
{
  fprintf (stderr, "\nUsage: mpirun -n Nr_of_nodes ./facs [options]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "\t-r reference Bloom filter to query against\n");
  fprintf (stderr, "\t-q FASTA/FASTQ file containing the query\n");
  fprintf (stderr, "\t-t threshold value\n");
  fprintf (stderr, "\t-f report output format, valid values are:\
           'json' and 'tsv'\n");
  fprintf (stderr, "\t-s sampling rate, default is 1 so it reads the whole\
           query file\n");
  fprintf (stderr, "\n");
/*finalize mpi process before each exit*/
  MPI_Finalize ();
  exit(1);
}
/*-------------------------------------*/
main (int argc, char **argv)
{
/*get MPI thead number and thread IDs*/
  int proc_num = 0, total_proc = 0;
/*----------MPI initialize------------*/
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &proc_num);
  MPI_Comm_size (MPI_COMM_WORLD, &total_proc);
/*------------variables----------------*/
  double tole_rate = 0, sampling_rate = 1;
  char *map_chunk = NULL, *bloom_filter = NULL, *target_path = NULL, *position = NULL, *position2 = NULL, *query = NULL, *report_fmt = "json";
  int opt=0;
  BIGCAST share=0, offset=0;
  char type = '@';
  gzFile zip = NULL;
  bloom *bl_2 = NEW (bloom);
  Queue *head = NEW (Queue), *tail = NEW (Queue), *head2 = head;
  head->location = NULL;
  head->next = tail;
/*------------get opt------------------*/
  
  if (argc<3)
  {
	if (proc_num == 0)
	{
        	return mpicheck_usage();
	}
  }
  while ((opt = getopt (argc, argv, "s:t:r:o:q:f:h")) != -1)
  {
      switch (opt)
      {
        case 't':
          tole_rate = atof(optarg);
          break;
        case 's':
          sampling_rate = atof(optarg);
          break;
        case 'o':
          target_path = optarg;
          break;
        case 'q':
          query = optarg;
          break;
        case 'r':
          bloom_filter = optarg;
          break;
        case 'f': // "json", "tsv" or none
          (optarg) && (report_fmt = optarg, 1);
          break;
        case 'h':
	  if (proc_num==0)
          	return mpicheck_usage();
          break;
        case '?':
          printf ("Unknown option: -%c\n", (char) optopt);
          if (proc_num==0)
	  	return mpicheck_usage();
          break;
      }
  }
  if (!bloom_filter && !query)
  {
	if (proc_num==0)
  		fprintf (stderr, "\nPlease, at least specify a bloom filter (-r) and a query file (-q)\n");
	MPI_Finalize ();
	exit (EXIT_FAILURE);
  }
  if (target_path == NULL)
  {
        target_path = argv[0];
  }
  if ((zip = gzopen (query, "rb")) < 0)
  {
	if (proc_num==0)
	{
        	fprintf(stderr, "%s\n", strerror(errno));
	}
		MPI_Finalize ();
        	exit(-1);
  }
  /*file type identification*/
  if (strstr (query, ".fastq") != NULL || strstr (query, ".fq") != NULL)
        type = '@';
  else
        type = '>';
  /*initialize emtpy string for query*/
  position = (char *) calloc ((ONEG + 1), sizeof (char));
  position2 = position;
  /*make structor for omp parallization*/
  F_set *File_head = make_list (bloom_filter, NULL);
  /*necessary initialization for python interface*/
  File_head->reads_num = 0;
  File_head->reads_contam = 0;
  File_head->hits = 0;
  File_head->all_k = 0;
  File_head->filename = bloom_filter;
  static char timestamp[40] = {0};
  // Get current timestamp, for benchmarking purposes (begin of run timestamp)
  isodate(timestamp);
  /*load a bloom filter*/
  load_bloom (File_head->filename, bl_2);
  /*load a default match cutoff according to k-mer*/
  if (tole_rate == 0)
  {
  	tole_rate = mco_suggestion (bl_2->k_mer);
  }
  /*pagelize chunk*/
  int page = getpagesize(); 
  int chunk = ONEG/page;
  /*get task share for each mpi thread*/
  share = struc_init (query,proc_num,total_proc,page);
  /*distribute task*/
  offset += share * proc_num;
  //printf ("share->%lld\n",share);
  //printf ("offset->%lld\n",offset);
  //printf ("chunk->%d\n",chunk);
  /*start processing*/
  while (share>0)
  {
	/*last piece must no bigger than chunk size*/
	if (share<=chunk)
	{
		chunk = share;
	}
	//printf ("share->%lld\n",share);
	//printf ("offset->%lld\n",offset);
	/*copy data from mapping area since standard read_process*/
	map_chunk = ammaping(query, offset, chunk, page);
	
	memcpy(position,map_chunk,chunk*page);
	/*remove fragments*/
	
	if (offset!=0)
	{
		position = remove_frag_head(position,type);
	}
	
	//if (proc_num!=total_proc || chunk<=ONEG)
	remove_frag_tail(position,type);
	//printf("reads num->%d\n",count_characters(position, '\n')/4);
	//printf ("real length->%d\n",strlen(position2));
	head = head2;
	head->next = tail;
	get_parainfo (position, head, type);

#pragma omp parallel
  	{
#pragma omp single nowait
	{
	  	while (head != tail)
		{
#pragma omp task firstprivate(head)
		{
	      		if (head->location != NULL)
                	{
				//printf("head->%0.25s\n",position);
				//printf ("position->%0.25s\n",head->location);
                        	read_process (bl_2, head, tail, File_head, sampling_rate, tole_rate, 'c', type);
			}
		}
	      	head = head->next;
		}
      	}	
	}
	offset+=chunk;
  	share-=chunk;
	position = position2;
  	munmap (map_chunk,chunk*page);
	if (share<chunk)
	{
     		memset (position,0,chunk*page);
	}
  }
  printf ("finish processing...\n");
  MPI_Barrier (MPI_COMM_WORLD);	//wait until all nodes finish
  gather (File_head, total_proc, proc_num);			//gather info from all nodes
  if (proc_num == 0)		
  {
  	char *result =  report(File_head, query, report_fmt, target_path, timestamp, prob_suggestion(bl_2->k_mer));
  	printf("%s\n",result);
  }
  MPI_Finalize ();
  return 0;
}
/*gather info from MPI thread*/
int gather (F_set *File_head, int total_proc, int proc_num)
{
  printf ("gathering...\n");
  BIGCAST temp, temp2, temp3, temp4;
  if (proc_num == 0)
  {
        // The master thread will need to receive all computations from all other threads.
  	MPI_Status status;
        // MPI_Recv(void *buf, int count, MPI_DAtatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
        // We need to go and receive the data from all other threads.
        // The arbitrary tag we choose is 1, for now.
  	int i = 0;
     	for (i = 1; i < total_proc; i++)
      	{
		MPI_Recv (&temp, 1, MPI_LONG_LONG_INT, i, 1, MPI_COMM_WORLD, &status);
	  	MPI_Recv (&temp2, 2, MPI_LONG_LONG_INT, i, 1, MPI_COMM_WORLD, &status);
	  	MPI_Recv (&temp3, 3, MPI_LONG_LONG_INT, i, 1, MPI_COMM_WORLD, &status);
		MPI_Recv (&temp4, 4, MPI_LONG_LONG_INT, i, 1, MPI_COMM_WORLD, &status);
	  	printf ("RECEIVED %lld from thread %d\n", temp, i);
	  	File_head->reads_num += temp;
	  	File_head->reads_contam += temp2;
	  	File_head->hits += temp3;
		File_head->all_k+=temp4;
		
      	}
  }
  else
  {
	temp = File_head->reads_num;
	temp2 = File_head->reads_contam;
	temp3 = File_head->hits;
	temp4 = File_head->all_k;
      	// We are finished with the results in this thread, and need to send the data to thread 1.
      	// MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
      	// The destination is thread 0, and the arbitrary tag we choose for now is 1.
        MPI_Send (&temp, 1, MPI_LONG_LONG_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send (&temp2, 2, MPI_LONG_LONG_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send (&temp3, 3, MPI_LONG_LONG_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send (&temp4, 4, MPI_LONG_LONG_INT, 0, 1, MPI_COMM_WORLD);
  }
  return 1;
}
/*distribute data by page*/
BIGCAST struc_init (char *filename, int proc_num, int total_proc, int page)
{
  BIGCAST total_size = get_size (filename);
  BIGCAST total_piece = total_size/page;
  BIGCAST share = 0;
  share = total_piece / total_proc;
  if (total_piece%total_proc!=0 && proc_num==(total_proc-1))
  {
        share += (total_piece % total_proc);
	/*
	if (total_size%page!=0)
	{
		share++;
	}
	*/
  } 
  return share; 
} 
/*partial mapping scheme*/
char *ammaping (char* source, BIGCAST offset, BIGCAST chunk, int page)
{
  struct stat statbuf;

  int fd = 0;
  char *sm = NULL;

  if ((fd = open (source, O_RDONLY | O_LARGEFILE)) < 0) {
      fprintf (stderr, "%s: %s\n", source, strerror (errno));
      MPI_Finalize ();
      exit (-1);
  }
  if (fstat (fd, &statbuf) < 0) {
      fprintf (stderr, "%s: %s\n", source, strerror (errno));
      MPI_Finalize ();
      exit (-1);
  } else if (statbuf.st_size == 0) {
      fprintf (stderr, "%s: %s\n", source, "File is empty");
      MPI_Finalize ();
      exit (-1);
  }
  /*offset by page*/
  sm = mmap (0, chunk*page, PROT_READ, MAP_SHARED | MAP_NORESERVE, fd, offset*page);                       

  if (MAP_FAILED == sm)
    {
      fprintf (stderr, "%s: %s\n", source, "Mapping error");
      MPI_Finalize ();
      exit (-1);
    }
  close (fd);
  return sm;
}
/*remove fragment in head area*/
char *remove_frag_head(char *position, char type)
{
	if (type == '@')
	{
		position = strstr (position,"\n+");
		position = strchr (strchr(position+1,'\n')+1,'\n')+1;
	}
	else
	{
		position = strchr (position, '>');
	}
return position;
}
/*remove fragment in tail area*/
void remove_frag_tail(char *position, char type)
{
	char *end = NULL;
	if (type == '@')
	{
		end = strrstr (position, "\n+");
        	end = bac_2_n (end - 1);
	}
	else
	{
		end = strrchr (end, '>') - 1;	
	}
	memset(end,0,strlen(end));
}

int count_characters(const char *str, char character)
{
    const char *p = str;
    int count = 0;

    do {
        if (*p == character)
            count++;
    } while (*(p++));

    return count;
}
