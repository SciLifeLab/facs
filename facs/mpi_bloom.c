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
#include "bloom.h"
#include "hashes.h"

#include<omp.h>
#include<mpi.h>

Queue *head, *head2, *tail;
void struc_init ();
/*-------------------------------------*/
int gather ();
/*-------------------------------------*/
char *ammaping (char *source);
/*-------------------------------------*/
mpi_main (int argc, char **argv)
{
  double tole_rate = 0, sampling_rate = 1;
  char *ref = NULL, *list = NULL, *target_path = NULL, *source = NULL, *report_fmt = "json";
  int opt, ntask = 0, mytask = 0, page = getpagesize (), total_piece;
  BIGCAST buffer=0, share=0, offset=0, chunk=0, total_size=0;
  bloom *bl_2 = NEW (bloom);
  Queue *head = NEW (Queue), *tail = NEW (Queue), *head2 = head;
  head->location = NULL;
  head->next = tail;
  while ((opt = getopt (argc, argv, "s:t:r:o:q:l:f:h")) != -1)
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
          source = optarg;
          break;
        case 'r':
          ref = optarg;
          break;
        case 'l':
          list = optarg;
          break;
        case 'f': // "json", "tsv" or none
          (optarg) && (report_fmt = optarg, 1);
          break;
        case 'h':
          return query_usage();
          break;
        case '?':
          printf ("Unknown option: -%c\n", (char) optopt);
          return query_usage();
          break;
      }
  }
  if (!target_path && !source)
  {
                fprintf (stderr, "\nPlease, at least specify a bloom filter (-r) and a query file (-q)\n");
                exit (-1);
  }
  if (target_path == NULL)
  {
                target_path = argv[0];
  } 
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &ntask);
  MPI_Comm_rank (MPI_COMM_WORLD, &mytask);
  struc_init ();		//structure init
  load_bloom (all_ref, bl_2);
  while (share > 0)
  {
  	position = ammaping (source);
  	get_parainfo (position);
   	head = head->next;
#pragma omp parallel
  	{
#pragma omp single nowait
	{
	  	while (head != tail)
		{
#pragma omp task firstprivate(head)
		{
			printf ("position->%0.10s\n", head->location);
	      		if (head->location != NULL)
                	{
                        	read_process (bl_2, head, tail, File_head, sampling_rate, tole_rate, mode, type);
			}
		}
	      	head = head->next;
		}
      	}	
	}
     		munmap (position, buffer * PAGE);
      		share -= buffer;
      		offset += buffer;
  }
  printf ("finish processing...\n");
  MPI_Barrier (MPI_COMM_WORLD);	//wait until all nodes finish
  gather ();			//gather info from all nodes
  if (mytask == 0)		
  {
	return report(File_head, query, report_fmt, target_path);
  }
  MPI_Finalize ();
  return 0;
}
/*-------------------------------------*/
void
struc_init (int *chunk, int *total_piece, BIGCAST *offset, BIGCAST *share, int page, int ntask, int mytask)
{
  get_size (source);		//get total size of file
  share = total_piece / ntask;	//every task gets an euqal piece
  if (total_piece % ntask != 0 && mytask == (ntask - 1))
  	share += (total_piece % ntask);	//last node tasks extra job
  offset = share * mytask;	//distribute the task
}
/*-------------------------------------*/
char *
ammaping (char *source)
{
  int src;
  char *sm;

  if ((src = open (source, O_RDONLY | O_LARGEFILE)) < 0)
    {
      perror (" open source ");
      exit (EXIT_FAILURE);
    }

  if (fstat (src, &statbuf) < 0)
    {
      perror (" fstat source ");
      exit (EXIT_FAILURE);
    }

  printf ("share->%d PAGES per node\n", share);

  if (share >= CHUNK)
    buffer = CHUNK;
  else
    buffer = share;
  printf ("total pieces->%d\n", total_piece);
  printf ("PAGE->%d\n", PAGE);
  printf ("node %d chunk size %d buffer size %d offset %d\n", mytask, CHUNK,
	  buffer, offset);

  sm = mmap (0, buffer * PAGE, PROT_READ, MAP_SHARED | MAP_NORESERVE, src, offset * PAGE);	//everytime we process a chunk of data

  //sm = mmap (0,share*PAGE, PROT_READ, MAP_SHARED | MAP_NORESERVE,src, offsetmytask*share*PAGE); //last time we process the rest

  if (MAP_FAILED == sm)
    {
      perror (" mmap source ");
      exit (EXIT_FAILURE);
    }

  return sm;
}

/*-------------------------------------*/
void
get_size (char *strFileName)
{
  stat (strFileName, &statbuf);
  PAGE = getpagesize ();	//get memory PAGE definition 
  total_piece = statbuf.st_size / PAGE;
  total_size = statbuf.st_size;
  CHUNK = 1000 * 1000 * 1000 * 1 / PAGE;	//1GB


  //if (statbuf.st_size % PAGE != 0)    //need one more page if total data is not a time number of a memory PAGE
  //extra_piece = statbuf.st_size % PAGE;
  //printf ("extra_piece->%d\n", extra_piece);
}

/*-------------------------------------*/
int
gather ()
{
  printf ("gathering...\n");
  if (mytask == 0)
    {
      // The master thread will need to receive all computations from all other threads.
      MPI_Status status;
      // MPI_Recv(void *buf, int count, MPI_DAtatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
      // We need to go and receive the data from all other threads.
      // The arbitrary tag we choose is 1, for now.
      int i = 0;
      for (i = 1; i < ntask; i++)
	{
	  BIGCAST temp, temp2, temp3;
	  MPI_Recv (&temp, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
	  MPI_Recv (&temp2, 5, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
	  MPI_Recv (&temp3, 7, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
	  printf ("RECEIVED %lld from thread %d\n", temp, i);
	  reads_num += temp;
	  reads_contam += temp2;
	  checky += temp3;
	}
    }
  else
    {
      // We are finished with the results in this thread, and need to send the data to thread 1.
      // MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
      // The destination is thread 0, and the arbitrary tag we choose for now is 1.
      MPI_Send (&reads_num, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send (&reads_contam, 5, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send (&checky, 7, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }
  return 1;
}

/*-------------------------------------*/

