#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <math.h>
#include <time.h>
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
#include "bloom.h"
#include "hashes.h"
/*-------------------------------------*/
//openMP library
#include<omp.h>
//MPICH/OPENMPI
#include<mpi.h>
/*-------------------------------------*/
//#define PERMS 0600
//#define NEW(type) (type *) malloc(sizeof(type))
/*-------------------------------------*/
int ntask = 0, mytask = 0;
/*-------------------------------------*/
long long total_piece, PAGE, buffer, share, offset, reads_num =
  0, reads_contam = 0, checky = 0, CHUNK, total_size = 0;
/*-------------------------------------*/
float error_rate, sampling_rate, contamination_rate, tole_rate;
/*-------------------------------------*/
int k_mer = 21, mode, mytask, ntask, type = 2, excel1, excel2, last_piece =
  0, extra_piece = 0;
int last = 0;
/*-------------------------------------*/
char *source, *all_ref, *position, *prefix;
/*-------------------------------------*/
Queue *head, *head2, *tail;
/*-------------------------------------*/
bloom *bl_2;
/*-------------------------------------*/
struct stat statbuf;
/*-------------------------------------*/
void list_init ();
void struc_init ();
void get_parainfo (char *full);
void get_size (char *strFileName);
void init (int argc, char **argv);
void fasta_process (bloom * bl, Queue * info);
void fastq_process (bloom * bl, Queue * info);
void evaluate (char *detail, char *filename);
void statistic_save (char *detail, char *filename);
/*-------------------------------------*/
int gather ();
int fastq_full_check (bloom * bl, char *p, int distance);
int fasta_full_check (bloom * bl, char *begin, char *next, char *model);
int fastq_read_check (char *begin, int length, char *model, bloom * bl);
int fasta_read_check (char *begin, char *next, char *model, bloom * bl);
/*-------------------------------------*/
char *jump (char *target);
char *ammaping (char *source);
//char* reallocate(Queue *info);
/*-------------------------------------*/
main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);

  MPI_Comm_size (MPI_COMM_WORLD, &ntask);

  MPI_Comm_rank (MPI_COMM_WORLD, &mytask);

  long sec, usec, i;

  struct timezone tz;

  struct timeval tv, tv2;

  if (mytask == 0)		//starting time
    {
      gettimeofday (&tv, &tz);
    }

  init (argc, argv);		//initialize 

  if (mode != 1 && mode != 2)
    {
      perror ("Mode select error.");
      return -1;
    }

  char *detail = (char *) malloc (1000 * 1000 * sizeof (char));

  memset (detail, 0, 1000 * 1000);

  strcat (detail, "query-->");
  strcat (detail, source);

  struc_init ();		//structure init

  load_bloom (all_ref, bl_2);

  k_mer = bl_2->k_mer;

  while (share > 0)
    {
      position = ammaping (source);

      list_init ();

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

		if (type == 1)
		  fasta_process (bl_2, head);
		else
		  fastq_process (bl_2, head);

	      }
	      head = head->next;

	    }
	}
      }
      munmap (position, buffer * PAGE);

      share -= buffer;

      offset += buffer;

      //head = head2;
    }
  printf ("finish processing...\n");

  MPI_Barrier (MPI_COMM_WORLD);	//wait until all threads finish jobs

  gather ();			//gather all matched and missed info

  if (mytask == 0)		//finishing time
    {
      gettimeofday (&tv2, &tz);

      sec = tv2.tv_sec - tv.tv_sec;

      usec = tv2.tv_usec - tv.tv_usec;

      printf ("total=%ld sec\n", sec);

      printf ("all->%d\n", excel1);

      printf ("excecuted->%d\n", excel2);

      evaluate (detail, all_ref);

      statistic_save (detail, source);
    }

  MPI_Finalize ();

  return 0;
}

/*-------------------------------------*/
void
init (int argc, char **argv)
{
  if (argc == 1 || !strcmp (argv[1], "-h") || !strcmp (argv[1], "-help"))
    {
      help ();
      check_help ();
    }
/*-------default-------*/
  mode = 1;
  k_mer = 21;
  tole_rate = 0.8;
  error_rate = 0.0005;
  sampling_rate = 1;
  prefix = NULL;

  int x;
  while ((x = getopt (argc, argv, "e:k:m:t:o:r:q:s:")) != -1)
    {
      //printf("optind: %d\n", optind);
      switch (x)
	{
	case 'e':
	  printf ("Error rate: \nThe argument of -e is %s\n", optarg);
	  (optarg) && ((error_rate = atof (optarg)), 1);
	  break;
	case 'k':
	  printf ("K_mer size: \nThe argument of -k is %s\n", optarg);
	  (optarg) && ((k_mer = atoi (optarg)), 1);
	  break;
	case 'm':
	  printf ("Mode : \nThe argument of -m is %s\n", optarg);
	  (optarg) && ((mode = atoi (optarg)), 1);
	  break;
	case 't':
	  printf ("Tolerant rate: \nThe argument of -t is %s\n", optarg);
	  (optarg) && ((tole_rate = atof (optarg)), 1);
	  break;
	case 's':
	  printf ("Sampling rate: \nThe argument of -s is %s\n", optarg);
	  (optarg) && ((sampling_rate = atof (optarg)), 1);
	  break;
	case 'o':
	  printf ("Output : \nThe argument of -o is %s\n", optarg);
	  (optarg) && ((prefix = optarg), 1);
	  break;
	case 'r':
	  printf ("Bloom list : \nThe argument of -r is %s\n", optarg);
	  (optarg) && (all_ref = optarg, 1);
	  break;
	case 'q':
	  printf ("Query : \nThe argument of -q is %s\n", optarg);
	  (optarg) && (source = optarg, 1);
	  break;
	case '?':
	  printf ("Unknown option: -%c\n", (char) optopt);
	  exit (0);
	}

    }

  if (strstr (source, ".fasta") || strstr (source, ".fna"))
    type = 1;

  if ((!all_ref) || (!source))
    exit (0);

  if (mode != 1 && mode != 2)
    {
      perror ("Mode select error.");
      exit (0);
    }

}

/*-------------------------------------*/
void
struc_init ()
{

  bl_2 = NEW (bloom);

  //head2 = head;

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
void
get_parainfo (char *full)
{
  printf ("distributing...\n");

  char *temp = full;

  int cores = omp_get_num_procs ();

  int offsett = buffer * PAGE / cores;

  last_piece = buffer * PAGE - (cores - 1) * offsett;

  int add = 0;

  printf ("last piece->%d\n", last_piece);

  Queue *pos = head;

  if (type == 1)
    {
      for (add = 0; add < cores; add++)
	{
	  Queue *x = NEW (Queue);
	  temp = strchr (full, '>');	//drop the possible fragment
	  if (add != 0)
	    temp = strchr (full + offsett * add, '>');
	  x->location = temp;
	  x->number = add;
	  x->next = pos->next;
	  pos->next = x;
	  pos = pos->next;
	}
    }				// end if

  else
    {
      for (add = 0; add < cores; add++)
	{
	  Queue *x = NEW (Queue);
	  if (add == 0 && *full != '@')
	    temp = strstr (full, "\n@") + 1;	//drop the fragment
	  if (add != 0)
	    temp = strstr (full + offsett * add, "\n@") + 1;
	  x->location = temp;
	  x->number = add;
	  x->next = pos->next;
	  pos->next = x;
	  pos = pos->next;
	}			//end else  
    }

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
	  long long temp, temp2, temp3;
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
void
fastq_process (bloom * bl, Queue * info)
{
  printf ("fastq processing...\n");
  char *p = info->location;
  char *next, *temp, *temp_piece = NULL;

  if (info->next != tail)
    next = info->next->location;

  else
    {
      printf ("last_piece  %d\n", last_piece);
      temp_piece = (char *) malloc ((last_piece + 1) * sizeof (char));
      memset (temp_piece, 0, last_piece + 1);
      memcpy (temp_piece, info->location, last_piece - PAGE);
      temp_piece[last_piece] = '\0';

      temp = temp_piece;
      while ((temp = strstr (temp, "\n@")))
	{
	  next = temp + 1;
	  temp++;
	}
      p = temp_piece;
    }

  while (p != next)
    {

      temp = jump (p);		//generate random number and judge if need to scan this read
      if (p != temp)
	{
	  p = temp;
	  continue;
	}
      if (p == '\0' || p == NULL)
	break;
#pragma omp atomic
      reads_num++;
      p = strchr (p, '\n') + 1;
      int distance = strchr (p, '\n') - p;

      if (!fastq_read_check (p, distance, "normal", bl))
#pragma omp atomic
	reads_contam++;
//printf("p->%s\n",p);
      p = strchr (p, '\n') + 1;
      p = strchr (p, '\n') + 1;
      p = strchr (p, '\n') + 1;
    }				// outside while
  printf ("finish process...\n");
  if (temp_piece)
    free (temp_piece);
}

/*-------------------------------------*/
int
fastq_read_check (char *begin, int length, char *model, bloom * bl)
{

//printf("fastq read check...\n");

  char *p = begin;
  int distance = length;
  int signal = 0, pre_kmer = 10, count_mis = 0;
  char *previous, *key = (char *) malloc (k_mer * sizeof (char) + 1);
  // int x = bloom_check(bl,"aaaaaaaaaaccccc");
  // printf("xxxxx->%d\n",x);
  while (distance > 0)
    {
      //printf("distance->%d\n",distance);
      if (signal == 1)
	break;

      if (distance >= k_mer)
	{
	  memcpy (key, p, sizeof (char) * k_mer);	//need to be tested
	  key[k_mer] = '\0';
	  previous = p;
	  p += k_mer;
	  distance -= k_mer;
	}

      else
	{
	  memcpy (key, previous + distance, sizeof (char) * k_mer);
	  p += (k_mer - distance);
	  signal = 1;
	}

      if (model == "reverse")
	rev_trans (key);

      //printf("key->%s\n",key);          
      if (mode == 1)
	{
	  if (bloom_check (bl, key))
	    {
	      //printf("hit\n");
	      return fastq_full_check (bl, begin, length);
	      //return 0;
	    }
	  //else
	  //printf("unhit\n");
	}

      else
	{
	  if (!bloom_check (bl, key))
	    {
	      //printf("unhit\n");
	      return fastq_full_check (bl, begin, length);
	      //return 0;
	    }
	  //else
	  //printf("hit\n");
	}
    }				// inner while

  free (key);

  if (model == "normal")	//use recursion to check the sequence forward and backward
    return fastq_read_check (begin, length, "reverse", bl);
  else
    return 1;
}

/*-------------------------------------*/
int
fastq_full_check (bloom * bl, char *p, int distance)
{

//printf("fastq full check...\n");

  int signal = 0, pre_kmer = -1;

  int count_mis = 0, label_m = 0, label_mis = 0, count = 0, match_s = 0;

  char *previous, *key = (char *) malloc (k_mer * sizeof (char) + 1);

//printf("k_mer->%din",k_mer);
  int length = distance;
#pragma omp atomic
  checky++;

  while (distance >= k_mer)
    {
      memcpy (key, p, sizeof (char) * k_mer);
      key[k_mer] = '\0';
      previous = p;
      p += 1;
      if (bloom_check (bl, key))
	{
	  count++;
	  if (pre_kmer == 1)
	    {
	      label_m++;
	      if (count < (k_mer-1))
		match_s++;
	      else
		{
		  match_s += count;
		  count = 0;
		}
	    }
	  else
	    {
	      label_m += k_mer;
	      match_s += k_mer - 1;
	    }
	  pre_kmer = 1;
	  //printf("%d----%d\n",label_m,label_mis);
	}
      else
	{
	  count = 0;
	  pre_kmer = 0;
	}
      distance--;
    }				// end while
  free (key);

  label_mis = length - label_m;

  if (((float) match_s / (float) (length)) >= tole_rate)
    return 0;
  else
    return 1;
}

/*-------------------------------------*/
void
fasta_process (bloom * bl, Queue * info)
{
  printf ("fasta processing...\n");
  char *p = info->location;
  char *next, *temp, *temp_next, *temp_piece = NULL;

  if (info->next != tail)
    next = info->next->location;

  else
    {
      last = 1;

      printf ("last_piece %d\n", last_piece);

      temp_piece = (char *) malloc ((last_piece + 1) * sizeof (char));

      memset (temp_piece, 0, last_piece + 1);

      memcpy (temp_piece, info->location, last_piece - PAGE);


      temp_piece[last_piece] = '\0';

      next = strrchr (temp_piece, '>');

      printf ("temp_piece->%0.30s\n", temp_piece);
      printf ("next->%0.30s\n", next);
      printf ("length->%d\n", strlen (temp_piece));
      printf ("test->%d\n", next - temp_piece);
      p = temp_piece;
      //printf ("p->%0.20s\n", p);
    }

  while (p != next)
    {
/*      
      if (last == 1){
         //printf ("p->%0.20s\n",p);
         printf("loop\n");
      }
*/
      //printf("here\n");
      temp = jump (p);		//generate random number and judge if need to scan this read
      if (p != temp)
	{
	  p = temp;
	  continue;
	}
#pragma omp atomic
      reads_num++;
      //if (last == 1)
      //printf("before\n");

      temp_next = strchr (p + 1, '>');

      //if (last == 1)
      //if (temp_next)
      //printf("temp_next->%0.20s\n",temp_next);
      //else
      //printf("boom\n");

      if (!temp_next)
	temp_next = next;

      if (!fasta_read_check (p, temp_next, "normal", bl))
	{
#pragma omp atomic
	  reads_contam++;
	}

      p = temp_next;
      /*
         if (last == 1)
         {
         printf ("p->%0.20s\n", p);
         printf ("temp_next->%0.20s\n",temp_next);
         if (p == next)
         printf("ja\n");
         }
       */
    }
  if (temp_piece)
    free (temp_piece);
}

/*-------------------------------------*/
int
fasta_read_check (char *begin, char *next, char *model, bloom * bl)
{
  //if (last == 1)
  //printf("fasta read check...\n");

//begin = strchr(begin+1,'\n')+1;

//printf("read_check->%0.10s\n",begin);

  char *p = strchr (begin + 1, '\n') + 1;

  char *start = p;

  char *key = (char *) malloc ((k_mer + 1) * sizeof (char));

  char *pre_key = (char *) malloc ((k_mer + 1) * sizeof (char));

  int n = 0, m = 0, count_enter = 0, result = 0;

  //int label_m = 0, label_mis = 0;

  key[k_mer] = '\0';


  while (p != next)
    {

      //if (last == 1)
      //{
      //printf("p->%0.30s\n",p);
      //printf("next->%0.30s\n",next);
      /*
         printf("next->%0.30s\n",next);
         char *mov=p;
         while (mov != next)
         {
         printf("(((->%0.30s\n",mov);
         mov = strchr(mov,'>')+1;
         }
       */
      //}

      while (n < k_mer)
	{
	  if (p[m] == '>' || p[m] == '\0')
	    {
	      m--;
	      break;
	    }

	  if (p[m] != '\r' && p[m] != '\n')
	    key[n++] = p[m];
	  else
	    count_enter++;
	  m++;
	}			//inner while

      if (m == 0)
	break;

      if (strlen (key) == k_mer)
	memcpy (pre_key, key, sizeof (char) * (k_mer + 1));

      else
	{
	  char *temp_key = (char *) malloc (k_mer * sizeof (char));

	  memcpy (temp_key, pre_key + strlen (key), k_mer - strlen (key));
	  memcpy (temp_key + k_mer - strlen (key), key,
		  sizeof (char) * (strlen (key) + 1));
	  free (key);
	  key = temp_key;

	}

      p += m;

      n = 0;

      m = 0;

      //if (last == 1)
      //printf("key->%s\n",key);

      if (model == "reverse")
	rev_trans (key);

      if (mode == 1)
	{
	  //printf("in\n");
	  if (bloom_check (bl, key))
	    {
	      //printf("in\n");
	      return fasta_full_check (bl, begin, next, model);
	      //return 0;
	    }
	  //else
	}			//outside if

      else
	{
	  if (!bloom_check (bl, key))
	    {
	      //printf("unhit\n");
	      return fasta_full_check (bl, begin, next, model);
	      //return 0;
	    }
	  //else
	  //printf("hit\n");
	}			//outside else
      memset (key, 0, k_mer);
    }				//outside while

  free (pre_key);
  free (key);

  if (model == "normal")	//use recursion to check the sequence forward and backward
    return fasta_read_check (begin, next, "reverse", bl);
  else
    {
//printf("in\n");
/*
if (mode==1)
label_mis+=(next-start-count_enter+1);
else
label_m+=(next-start-count_enter+1);
*/
//printf("one read finish...\n");
      return 1;
    }
}

/*-------------------------------------*/
int
fasta_full_check (bloom * bl, char *begin, char *next, char *model)
{

#pragma omp atomic
  checky++;

  int label_m = 0, label_mis = 0, match_s = 0, count = 0;

  //printf ("fasta full check...\n");

  int n = 0, m = 0, count_enter = 0, pre_kmer = -1;

  char *key = (char *) malloc ((k_mer + 1) * sizeof (char));

  //printf("%0.10s\n",begin);

  begin = strchr (begin + 1, '\n') + 1;

  char *p = begin;

  while (p != next)
    {
      if (*p == '\n')
	count_enter++;
      p++;
    }

  p = begin;

  while (*p != '>' && *p != '\0')
    {
      while (n < k_mer)
	{
	  //printf("k_mer...\n");
	  if (p[m] == '>' || p[m] == '\0')
	    {
	      m--;
	      break;
	    }

	  if (p[m] != '\r' && p[m] != '\n')
	    key[n++] = p[m];

	  m++;
	}
      key[n] = '\0';

      if (model == "reverse")
	rev_trans (key);

      if (strlen (key) == k_mer)
	if (bloom_check (bl, key))
	  {
	    count++;
	    if (pre_kmer == 1)
	      {
		label_m++;
		if (count < (k_mer-1))
		  match_s++;
		else
		  {
		    match_s += count;
		    count = 0;
		  }
	      }
	    else
	      {
		label_m += k_mer;
		match_s += k_mer - 1;
	      }
	    pre_kmer = 1;
	    //printf("%d----%d\n",label_m,label_mis);
	  }
	else
	  {
	    count = 0;
	    pre_kmer = 0;
	  }

      p++;
      if (p[0] == '\n')
	p++;
      n = 0;
      m = 0;
    }				// end of while

  if (((float) (match_s) / (float) (next - begin - count_enter)) >= (tole_rate))	//match >tole_rate considered as contaminated
    return 0;
  else
    return 1;
}

/*-------------------------------------*/
void
evaluate (char *detail, char *filename)
{
  char buffer[100] = { 0 };
  printf ("all->%d\n", reads_num);
  printf ("contam->%d\n", reads_contam);
  printf ("possbile->%d\n", checky);
  contamination_rate = (double) (reads_contam) / (double) (reads_num);
  if ((mode == 1 && contamination_rate == 0)
      || (mode == 2 && contamination_rate == 1))
    printf ("clean data...\n");
  else if (mode == 1)
    printf ("contamination rate->%f\n", contamination_rate);
  else
    printf ("contamination rate->%f\n", 1 - contamination_rate);

  strcat (detail, "\nxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  strcat (detail, "bloom->");
  strcat (detail, filename);
  strcat (detail, "   \n");
  sprintf (buffer, "all->%d\n", reads_num);
  strcat (detail, buffer);
  memset (buffer, 0, 100);
  sprintf (buffer, "contam->%d\n", reads_contam);
  strcat (detail, buffer);
  memset (buffer, 0, 100);
  sprintf (buffer, "possbile->%d\n", checky);
  strcat (detail, buffer);
  memset (buffer, 0, 100);
  sprintf (buffer, "contamination rate->%f", contamination_rate);
  strcat (detail, buffer);
  memset (buffer, 0, 100);
  reads_num = 0;
  reads_contam = 0;
  contamination_rate = 0;
}

/*-------------------------------------*/
char *
jump (char *target)
{
  //printf("here\n");
  float seed = rand () % 10;

  if (seed >= (float) sampling_rate * 10)
    {

      char *point;

      if (type == 1)
	point = strchr (target + 1, '>');	//point to >
      else
	point = strstr (target + 1, "\n@") + 1;	//point to @

      if (point)
	target = point;
    }
  return target;
}

/*-------------------------------------*/
void
statistic_save (char *detail, char *filename)
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

  printf ("bloom name->%s\n", save_file);

  write_result (save_file, detail);
}

void
list_init ()
{
  head = NEW (Queue);

  tail = NEW (Queue);

  head->next = tail;
}

/*

char* reallocate(Queue *info)
{
  char *temp_piece = (char *) malloc (last_piece*sizeof(char)+1);
  memcpy (temp_piece,info->location,last_piece);
  //rintf("here\n");
  temp_piece[last_piece]='\0';
  return temp_piece;
}

*/
