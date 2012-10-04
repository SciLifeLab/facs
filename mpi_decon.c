#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <stdio.h>
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
#include<mpi.h>
/*-------------------------------------*/
#define PERMS 0600
#define NEW(type) (type *) malloc(sizeof(type))
/*-------------------------------------*/
long sec, usec, i;
struct timezone tz;
struct timeval tv, tv2;
/*-------------------------------------*/
int ntask = 0, mytask = 0;
/*-------------------------------------*/
long long total_piece, PAGE, buffer, share, offset, CHUNK;
/*-------------------------------------*/
float error_rate, tole_rate, contamination_rate;
/*-------------------------------------*/
int k_mer, mode, mytask, ntask, type = 0, count1 = 0, last_piece;
/*-------------------------------------*/
char *source, *all_ref, *position, *prefix, *clean, *contam, *clean2,
  *contam2, *merge_contam, *merge_clean;
/*-------------------------------------*/
Queue *head, *head2, *tail;
/*-------------------------------------*/
bloom *bl_2;
/*-------------------------------------*/
//F_set *File_head;
/*-------------------------------------*/
void list_init ();
void task_init ();
void struc_init ();
void distribut_init ();
void get_parainfo (char *full);
void save_result (char *source);
void get_size (char *strFileName);
void init (int number, char **orders);
void para_init (int argc, char **argv);
void fasta_process (bloom * bl, Queue * info);
void fastq_process (bloom * bl, Queue * info);
void cat (char *clean, char *contam, char *source);
/*-------------------------------------*/
//int gather();
int fastq_full_check (bloom * bl, char *p, int distance);
int fasta_full_check (bloom * bl, char *begin, char *next, char *model);
int fastq_read_check (char *begin, int length, char *model, bloom * bl);
int fasta_read_check (char *begin, char *next, char *model, bloom * bl);
/*-------------------------------------*/
char *ammaping (char *source);
/*-------------------------------------*/
main (int argc, char **argv)
{

  para_init (argc, argv);

  init (argc, argv);		//initialize 

  distribut_init ();

  struc_init ();		//structure init

  task_init ();

  load_bloom (all_ref, bl_2);

  k_mer = bl_2->k_mer;

  while (share > 0)
    {

      if (strstr (source, ".fifo"))
	position = large_load (source);
      else
	position = ammaping (source);
      list_init ();
      get_parainfo (position);

      printf ("xxxxxxxxxxxxxaxxxxxxxxxxxxxxxxxx\n");

#pragma omp parallel
      {
#pragma omp single nowait
	{
	  while (head != tail)
	    {
#pragma omp task firstprivate(head)
	      {
		if (head->location)
		  {
		    printf ("location->%0.20s\n", head->location);

		    if (type == 1)
		      fasta_process (bl_2, head);
		    else
		      fastq_process (bl_2, head);

		  }
	      }
	      head = head->next;
	    }
	}			// End of single - no implied barrier (nowait)
      }				// End of parallel region - implied barrier

      munmap (position, buffer * PAGE);

      share -= buffer;

      offset += buffer;

      save_result (source);	//save result and memset clean/contam

      count1++;

    }				//end while

  MPI_Barrier (MPI_COMM_WORLD);

  printf ("finish processing...\n");

  cat (merge_clean, merge_contam, source);

  if (mytask == 0)		//finishing time
    {

      gettimeofday (&tv2, &tz);

      sec = tv2.tv_sec - tv.tv_sec;

      usec = tv2.tv_usec - tv.tv_usec;

      printf ("total=%ld sec\n", sec);

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
  prefix = NULL;

  int x;
  while ((x = getopt (argc, argv, "e:k:m:t:o:r:q:")) != -1)
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

  clean = (char *) malloc (CHUNK * PAGE * sizeof (char));

  contam = (char *) malloc (CHUNK * PAGE * sizeof (char));

  merge_contam = (char *) malloc (1000 * sizeof (char));

  merge_clean = (char *) malloc (1000 * sizeof (char));

  clean2 = clean;

  contam2 = contam;

  memset (merge_contam, 0, 1000);
  memset (merge_clean, 0, 1000);

  printf ("lo->%d\n", CHUNK * PAGE);
}

/*-------------------------------------*/
char *
ammaping (char *source)
{
  int src;
  char *sm;
  struct stat statbuf;

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

  printf ("task->%d\n", mytask);

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
  struct stat statbuf;
  stat (strFileName, &statbuf);
  PAGE = getpagesize ();	//get memory PAGE definition 
  total_piece = statbuf.st_size / PAGE;
  CHUNK = 1000 * 1000 * 1000 * 1 / PAGE;	//1GB
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

  Queue *pos = head;

  if (type == 0)
    if (*full == '>')
      type = 1;
    else if (*full == '@')
      type = 2;

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
  head = head->next;
}

/*-------------------------------------*/
void
fastq_process (bloom * bl, Queue * info)
{
  printf ("fastq processing...\n");

  char *p = info->location;
  char *next, *temp, *temp_start, *temp_end, *temp_piece = NULL;

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

      temp_start = p;

      if (p == '\0' || p == NULL)
	break;

      p = strchr (p, '\n') + 1;

      temp_end = strstr (p, "\n@");

      int result = fastq_read_check (p, strchr (p, '\n') - p, "normal", bl);

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
int
fastq_read_check (char *begin, int length, char *model, bloom * bl)
{

//printf("fastq read check...\n");

  char *p = begin;
  int distance = length;
  int signal = 0, pre_kmer = 10, count_mis = 0;
  char *previous, *key = (char *) malloc (k_mer * sizeof (char) + 1);

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

      if (mode == 1)
	{
	  if (bloom_check (bl, key))
	    {
	      return fastq_full_check (bl, begin, length);
	    }
	}

      else
	{
	  if (!bloom_check (bl, key))
	    {
	      return fastq_full_check (bl, begin, length);
	    }
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

  int length = distance;

  int signal = 0, pre_kmer = -1;

  int count_mis = 0, label_m = 0, label_mis = 0, count = 0, match_s = 0;

  char *previous, *key = (char *) malloc (k_mer * sizeof (char) + 1);

  while (distance >= k_mer)
    {
      memcpy (key, p, sizeof (char) * k_mer);
      key[k_mer] = '\0';
      previous = p;
      p += 1;
      //printf("key->%s\n",key);
      if (bloom_check (bl, key))
	{
	  count++;
	  if (pre_kmer == 1)
	    {
	      label_m++;
	      if (count < 20)
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

  if (label_m == 0 && label_mis == 0)
    return -1;
  if (((float) match_s / (float) length) < tole_rate)
    return 0;
  else
    return 1;
}

/*-------------------------------------*/
void
fasta_process (bloom * bl, Queue * info)
{
  printf ("fasta processing...\n");

  char *temp, *next, *p = info->location, *temp_piece = NULL;

  if (info->next != tail)
    next = info->next->location;
  else
    {
      temp_piece = (char *) malloc ((last_piece + 1) * sizeof (char));
      memcpy (temp_piece, info->location, last_piece - PAGE);
      temp_piece[last_piece] = '\0';
      next = strrchr (temp_piece, '>');
      p = temp_piece;
    }

  while (p != next)
    {
      temp = strchr (p + 1, '>');
      if (!temp)
	temp = next;
      int result = fasta_read_check (p, temp, "normal", bl);

      if (result)
	{
#pragma omp critical
	  {

	    memcpy (clean, p, temp - p);
	    clean += (temp - p);

	  }
	}
      else if (result == 0)
	{
#pragma omp critical
	  {
	    memcpy (contam, p, temp - p);
	    contam += (temp - p);
	  }
	}
      p = temp;
    }
  if (temp_piece)
    free (temp_piece);
}

/*-------------------------------------*/
int
fasta_read_check (char *begin, char *next, char *model, bloom * bl)
{

  char *p = strchr (begin + 1, '\n') + 1;

  if (!p || *p == '>')
    return 1;

  char *start = p;

  int n, m, label_m, label_mis, count_enter;

  char *key = (char *) malloc ((k_mer + 1) * sizeof (char));

  char *pre_key = (char *) malloc ((k_mer + 1) * sizeof (char));

  key[k_mer] = '\0';

  while (p != next)
    {
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

      if (model == "reverse")
	rev_trans (key);

      if (mode == 1)
	{
	  if (bloom_check (bl, key))
	    {
	      return fasta_full_check (bl, begin, next, model);
	    }
	}			//outside if

      else
	{
	  if (!bloom_check (bl, key))
	    {
	      return fasta_full_check (bl, begin, next, model);
	    }
	}			//outside else
      memset (key, 0, k_mer);
    }				//outside while

  if (model == "normal")	//use recursion to check the sequence forward and backward
    return fasta_read_check (begin, next, "reverse", bl);
  else
    return 1;
}

/*-------------------------------------*/
int
fasta_full_check (bloom * bl, char *begin, char *next, char *model)
{
  int label_m = 0, label_mis = 0, match_s = 0, count = 0;

  int n = 0, m = 0, count_enter = 0, pre_kmer = -1;

  char *key = (char *) malloc ((k_mer + 1) * sizeof (char));

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
		if (count < 20)
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
save_result (char *source)
{
  char *clean_name = (char *) malloc (200 * sizeof (char));
  char *contam_name = (char *) malloc (200 * sizeof (char));

  memset (clean_name, 0, 200);
  memset (contam_name, 0, 200);

  char *position, *position2;
  char *so = source;

  while (strchr (source, '/'))
    {
      position = strchr (source, '/');
      source = position + 1;
    }
  /*--series number--*/
  sprintf (clean_name, "%d", mytask * 10 + count1);
  sprintf (contam_name, "%d", mytask * 10 + count1);
  /*--tag--*/
  strcat (clean_name, "_clean");
  strcat (contam_name, "_contam");
  /*--tail2--*/
  if (type == 1)
    position = ".fasta";
  else
    position = ".fastq";
  strcat (clean_name, position);
  strcat (contam_name, position);

  //printf("clean name->%s\n",clean_name);
  //printf("contam name->%s\n",contam_name);

  *clean = '\0';
  *contam = '\0';

  write_result (clean_name, clean2);

  write_result (contam_name, contam2);

  strcat (merge_contam, "  ");

  strcat (merge_clean, "  ");

  strcat (merge_contam, clean_name);

  strcat (merge_clean, contam_name);

  //printf("merge->%s\n",merge_clean);

  //printf("merge->%s\n",merge_contam);

  free (clean_name);

  free (contam_name);

  memset (clean2, 0, CHUNK * PAGE * sizeof (char));
  memset (contam2, 0, CHUNK * PAGE * sizeof (char));

  clean = clean2;
  contam = contam2;
}

/*-------------------------------------*/
void
cat (char *clean, char *contam, char *source)
{

  printf ("clean->%s\n", clean);
  printf ("contam->%s\n", contam);

  if (mytask == 0)
    {
      char *so;
      ((so = strrchr (source, '/'))) && (so + 1) || (so = NULL);

      char *cat_clean = (char *) malloc (1000 * 1000 * sizeof (char));
      char *cat_contam = (char *) malloc (1000 * 1000 * sizeof (char));
      char *rm_clean = (char *) malloc (1000 * 1000 * sizeof (char));
      char *rm_contam = (char *) malloc (1000 * 1000 * sizeof (char));
      char *temp = (char *) malloc (1000 * sizeof (char));
      char *temp2 = (char *) malloc (1000 * sizeof (char));
      memset (cat_clean, 0, 1000 * 1000);
      memset (cat_contam, 0, 1000 * 1000);
      memset (rm_clean, 0, 1000 * 1000);
      memset (rm_contam, 0, 1000 * 1000);
      strcat (cat_clean, "cat ");
      strcat (cat_contam, "cat ");
      strcat (rm_clean, "rm ");
      strcat (rm_contam, "rm ");
      // The master thread will need to receive all computations from all other threads.
      MPI_Status status;
      // MPI_Recv(void *buf, int count, MPI_DAtatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
      // We need to go and receive the data from all other threads.
      // The arbitrary tag we choose is 1, for now.
      int i = 0;
      for (i = 1; i < ntask; i++)
	{
	  char *temp = (char *) malloc (1000 * sizeof (char)), *temp2 =
	    (char *) malloc (1000 * sizeof (char));
	  MPI_Recv (temp, 1000, MPI_CHAR, i, 1, MPI_COMM_WORLD, &status);
	  MPI_Recv (temp2, 1000, MPI_CHAR, i, 2, MPI_COMM_WORLD, &status);
	  printf ("RECEIVED %s from thread %d\n", temp, i);
	  printf ("RECEIVED %s from thread %d\n", temp2, i);
	  strcat (cat_clean, temp);
	  strcat (cat_contam, temp2);
	  strcat (rm_clean, temp);
	  strcat (rm_contam, temp2);
	  strcat (cat_clean, "  ");
	  strcat (cat_contam, "  ");
	  strcat (rm_clean, "  ");
	  strcat (rm_contam, "  ");
	  memset (temp, 0, 1000);
	  memset (temp2, 0, 1000);
	  printf ("clean order->%s\ncontam order->%s\n", cat_clean,
		  cat_contam);
	}
      strcat (cat_clean, clean);	//  for node 0
      strcat (cat_contam, contam);	//  for node 0
      strcat (rm_clean, clean);
      strcat (rm_contam, contam);

      strcat (cat_clean, " > ");
      strcat (cat_contam, " > ");

      if (prefix && so)
	{
	  strcat (cat_clean, prefix);
	  strcat (cat_contam, prefix);
	  strncat (cat_clean, so, strrchr (source, '.') - so);
	  strncat (cat_contam, so, strrchr (source, '.') - so);
	}
      else if (prefix && !so)
	{
	  strcat (cat_clean, prefix);
	  strcat (cat_contam, prefix);
	  strncat (cat_clean, source, strrchr (source, '.') - source);
	  strncat (cat_contam, source, strrchr (source, '.') - source);
	}
      else if (!prefix)
	{
	  strncat (cat_clean, source, strrchr (source, '.') - source);
	  strncat (cat_contam, source, strrchr (source, '.') - source);
	}

      strcat (cat_clean, "_clean_");
      strcat (cat_contam, "_contam_");

      if (type == 1)
	{
	  strcat (cat_clean, ".fasta");
	  strcat (cat_contam, ".fasta");
	}
      else
	{
	  strcat (cat_clean, ".fastq");
	  strcat (cat_contam, ".fastq");
	}
      printf ("final order-->\n%s\n%s\n", cat_clean, cat_contam);
      printf ("final remove-->\n%s\n%s\n", rm_clean, rm_contam);

      system (cat_clean);
      system (cat_contam);

      system (rm_clean);
      system (rm_contam);

    }
  else
    {
      MPI_Send (clean, strlen (clean), MPI_CHAR, 0, 1, MPI_COMM_WORLD);
      MPI_Send (contam, strlen (contam), MPI_CHAR, 0, 2, MPI_COMM_WORLD);
      printf ("sent...\n");
    }
}

/*-------------------------------------*/
void
para_init (int argc, char **argv)
{
  MPI_Init (&argc, &argv);

  MPI_Comm_size (MPI_COMM_WORLD, &ntask);

  MPI_Comm_rank (MPI_COMM_WORLD, &mytask);

  if (mytask == 0)		//starting time
    {
      gettimeofday (&tv, &tz);
    }
}

/*-------------------------------------*/
void
distribut_init ()
{
  get_size (source);		//get total size of file
  share = total_piece / ntask;	//every task gets an euqal piece
  if (total_piece % ntask != 0 && mytask == (ntask - 1))
    share += (total_piece % ntask);	//last node tasks extra job
  offset = share * mytask;
}

/*-------------------------------------*/
void
task_init ()
{
  bl_2 = NEW (bloom);
}

/*-------------------------------------*/
void
list_init ()
{
  head = NEW (Queue);
  tail = NEW (Queue);
  head->next = tail;
  head2 = head;
}
