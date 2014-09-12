#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64

#include <math.h>
#include <time.h>
#include <errno.h>
#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <string.h>

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>

#include "tool.h"
#include "prob.h"
#include "bloom.h"
#include "query.h"
#include "hashes.h"
#include "mpi_bloom.h"
#include <omp.h>
#include <mpi.h>

char *temp = NULL;
/*-------------------------------------*/
static int mpicheck_usage (void)
{
  fprintf (stderr, "\nUsage: mpirun -n Nr_of_nodes ./facs [options]\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "\t-r reference Bloom filter to query against\n");
  fprintf (stderr, "\t-q FASTA/FASTQ file containing the query\n");
  fprintf (stderr, "\t-l input list containing all Bloom filters,\
           one per line\n");
  fprintf (stderr, "\t-t threshold value\n");
  fprintf (stderr, "\t-f report output format, valid values are:\
           'json' and 'tsv'\n");
  fprintf (stderr, "\t-s sampling rate, default is 1 so it reads the whole\
           query file\n");
  fprintf (stderr, "\n");
  exit(1);
}
/*-------------------------------------*/
main (int argc, char **argv)
{
  int proc_num = 0, total_proc = 0;
/*----------MPI initialize------------*/
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &proc_num);
  MPI_Comm_size (MPI_COMM_WORLD, &total_proc);
/*------------variables----------------*/
  double tole_rate = 0, sampling_rate = 1;
  char *bloom_filter = NULL,  *target_path = NULL, *position = NULL, *query = NULL, *report_fmt = "json";
  int opt=0;
  BIGCAST share=0, offset=0;
  char type = '@';
  gzFile zip = NULL;
  static char timestamp[40] = {0};

  // Get current timestamp, for benchmarking purposes (begin of run timestamp)
  isodate(timestamp);
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
        	exit(EXIT_FAILURE);
  }
  if (strstr (query, ".fastq") != NULL || strstr (query, ".fq") != NULL)
        type = '@';
  else
        type = '>';
  /*initialize emtpy string for query*/
  position = (char *) calloc ((ONEG + 1), sizeof (char));
  share = struc_init (query,proc_num,total_proc);
  F_set *File_head = make_list (bloom_filter, NULL);
  File_head->reads_num = 0;
  File_head->reads_contam = 0;
  File_head->hits = 0;
  File_head->all_k = 0;
  File_head->filename = bloom_filter;
  load_bloom (File_head->filename, bl_2);	//load a bloom filter
  if (tole_rate == 0)
  {
  	tole_rate = mco_suggestion (bl_2->k_mer);
  }
  
  while (offset<share)
  {
        offset+= gz_mpi (zip, offset+share*proc_num, share-offset, position, type);

	head = head2;
	head->next = tail;
	if (temp)
	position = temp;
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
                        	read_process (bl_2, head, tail, File_head, sampling_rate, tole_rate, 'c', type);
			}
		}
	      	head = head->next;
		}
      	}	
	}

     		memset (position, 0, strlen(position));
  }
  
  printf ("finish processing...\n");
  MPI_Barrier (MPI_COMM_WORLD);	//wait until all nodes finish
  gather (File_head, total_proc, proc_num);			//gather info from all nodes
  if (proc_num == 0)		
  {
  	char *result =  report(File_head, query, report_fmt, target_path, timestamp, prob_suggestion(bl_2->k_mer), total_proc);
  	printf("%s\n",result);
  }
  MPI_Finalize ();
  return 0;
}
/*-------------------------------------*/
BIGCAST struc_init (char *filename, int proc_num, int total_proc)
{
  
  BIGCAST total_size = get_size(filename);
  BIGCAST share = 0;
  share = total_size / (BIGCAST)total_proc;	//every task gets an euqal piece
  if (total_size%total_proc!=0 && proc_num==(total_proc-1))
  {
  	share += (total_size % total_proc);	//last node takes extra job
  }
  return share;
}
/*-------------------------------------*/
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

/*-------------------------------------*/
BIGCAST gz_mpi (gzFile zip, BIGCAST offset, BIGCAST left, char *data, char type)
{
  char *start=NULL, *end = NULL;
  BIGCAST complete = 0;
  gzseek (zip, offset, SEEK_SET);
  if (left>(ONEG))
  {
  	gzread (zip, data, ONEG);
  	//complete = 2*ONEG;
  }
  else
  {
	gzread (zip, data, left);
  	complete = left;
  }
  if (type == '@')
  {
	if (offset!=0)
	{
		start = strstr (data,"\n+");
		start = strchr (strchr(start+1,'\n')+1,'\n')+1;
		temp = start;
        }
	end = strrstr (data, "\n+");
        end = move_start_point (end - 1);
  }
  else
  {
	if (offset!=0)
	{
		start = strchr (data, '>');
	}
        end = strrchr (data, '>') - 1;
  }
  if (start)
  {
	complete -= (start - data);
	//*data = start;
  }
  if (end)
  {
        complete += (end - data);
        memset (end, 0, strlen (end));
  }
  return complete;
}
/*-------------------------------------*/
