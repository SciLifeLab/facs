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
//#include<mpi.h>
/*-------------------------------------*/
#define PERMS 0600
#define NEW(type) (type *) malloc(sizeof(type))
/*-------------------------------------*/
long long share, offset, CHUNK, total_size;
/*-------------------------------------*/
float error_rate, tole_rate, contamination_rate;
/*-------------------------------------*/
int k_mer, mode, mytask, ntask, type = 2;
/*-------------------------------------*/
char *source, *all_ref, *position, *prefix, *clean, *contam, *clean2,
  *contam2;
/*-------------------------------------*/
Queue *head, *tail;
/*-------------------------------------*/
bloom *bl_2;
/*-------------------------------------*/
struct stat statbuf;
/*-------------------------------------*/
void struc_init ();
void get_parainfo (char *full);
void get_size (char *strFileName);
void init (int argc, char **argv);
void fasta_process (bloom * bl, Queue * info);
void fastq_process (bloom * bl, Queue * info);
void save_result (char *source, char *obj_file);
/*-------------------------------------*/
char *mmaping (char *source);
/*-------------------------------------*/
int fastq_full_check (bloom * bl, char *p, int distance);
int fasta_full_check (bloom * bl, char *begin, char *next, char *model);
int fastq_read_check (char *begin, int length, char *model, bloom * bl);
int fasta_read_check (char *begin, char *next, char *model, bloom * bl);
/*-------------------------------------*/
main (int argc, char **argv)
{

  long sec, usec, i;

  struct timezone tz;

  struct timeval tv, tv2;

  gettimeofday (&tv, &tz);

  Queue *head1;

  init (argc, argv);		//initialize 

  struc_init ();		//structure init   

  position = mmaping (source);

  get_parainfo (position);

  char *detail = (char *) malloc (1000 * 1000 * sizeof (char));

  memset (detail, 0, 1000 * 1000);

  load_bloom (all_ref, bl_2);

      k_mer = bl_2->k_mer;

#pragma omp parallel
      {
#pragma omp single nowait
	{
	  while (head != tail)
	    {
#pragma omp task firstprivate(head)
	      {
		if (head->location)
		  if (type == 1)
		    fasta_process (bl_2, head);
		  else
		    fastq_process (bl_2, head);
	      }
	      head = head->next;
	    }
	}			// End of single - no implied barrier (nowait)
      }				// End of parallel region - implied barrier

      save_result (source, all_ref);

  munmap (position, statbuf.st_size);

  printf ("finish processing...\n");

  gettimeofday (&tv2, &tz);

  sec = tv2.tv_sec - tv.tv_sec;

  usec = tv2.tv_usec - tv.tv_usec;

  printf ("total=%ld sec\n", sec);

  return 0;
}

/*-------------------------------------*/
void
init (int argc, char **argv)
{
  if (argc == 1 || !strcmp (argv[1], "-h") || !strcmp (argv[1], "-help"))
    {
      help ();
      remove_help ();
    }
/*-------default-------*/
  mode = 1;
  k_mer = 21;
  tole_rate = 0.8;
  error_rate = 0.0005;
  prefix = NULL;
/*-------default-------*/
  int x;
  while ((x = getopt (argc, argv, "e:k:m:t:o:r:q:")) != -1)
    {
      //printf("optind: %d\n", optind);
      switch (x)
	{
	case 'e':
	  //printf ("Error rate: \nThe argument of -e is %s\n", optarg);
	  (optarg) && ((error_rate = atof (optarg)), 1);
	  break;
	case 'k':
	  //printf ("K_mer size: \nThe argument of -k is %s\n", optarg);
	  (optarg) && ((k_mer = atoi (optarg)), 1);
	  break;
	case 'm':
	  //printf ("Mode : \nThe argument of -m is %s\n", optarg);
	  (optarg) && ((mode = atoi (optarg)), 1);
	  break;
	case 't':
	  //printf ("Tolerant rate: \nThe argument of -t is %s\n", optarg);
	  (optarg) && ((tole_rate = atof (optarg)), 1);
	  break;
	case 'o':
	  //printf ("Out : \nThe argument of -o is %s\n", optarg);
	  (optarg) && ((prefix = optarg), 1);
	  break;
	case 'r':
	  //printf ("Bloom list : \nThe argument of -r is %s\n", optarg);
	  (optarg) && ((all_ref = optarg), 1);
	  break;
	case 'q':
	  //printf ("Query : \nThe argument of -q is %s\n", optarg);
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

  printf ("all_ref->%s\n", all_ref);

}

/*-------------------------------------*/
void
struc_init ()
{

  bl_2 = NEW (bloom);

  head = NEW (Queue);

  tail = NEW (Queue);

  head->next = tail;

}

/*-------------------------------------*/
char *
mmaping (char *source)
{
  int src;
  char *sm;

  if ((src = open (source, O_RDONLY)) < 0)
    {
      perror (" open source ");
      exit (EXIT_FAILURE);
    }

  if (fstat (src, &statbuf) < 0)
    {
      perror (" fstat source ");
      exit (EXIT_FAILURE);
    }

  sm =
    mmap (0, (size_t) statbuf.st_size, PROT_READ, MAP_SHARED | MAP_NORESERVE,
	  src, 0);

  if (MAP_FAILED == sm)
    {
      perror (" mmap source ");
      exit (EXIT_FAILURE);
    }

  clean = (char *) malloc (statbuf.st_size * sizeof (char));

  contam = (char *) malloc (statbuf.st_size * sizeof (char));

  clean2 = clean;

  contam2 = contam;

  return sm;
}

/*-------------------------------------*/
void
get_parainfo (char *full)
{
  printf ("distributing...\n");

  char *temp=full;

  int cores = omp_get_num_procs ();

  int offsett = statbuf.st_size / cores;

  int add = 0;

  printf ("task->%d\n", offsett);

  Queue *pos = head;

  if (type == 1)
    {
      for (add = 0; add < cores; add++)
	{
	  Queue *x = NEW (Queue);

          if (add == 0 && *full != '>')

	  temp = strchr (full, '>');	//drop the possible fragment

	  if (add != 0)

	    temp = strchr (full + offsett, '>');

	  //printf ("full->%0.20s\n", full);

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

          //printf("offset->%d\n",offsett*add);

	  if (add != 0)

	    temp = strstr (full + offsett*add, "\n@");

          if (temp)
          temp++;

	  x->location = temp;

	  x->number = add;

	  x->next = pos->next;

	  pos->next = x;

	  pos = pos->next;
	}			//end else  

    }

  return;
}
/*-------------------------------------*/
void
fastq_process (bloom * bl, Queue * info)
{
  printf ("fastq processing...\n");

  int read_num = 0;
  char *p = info->location;
  char *next, *temp_start, *temp_end, *temp_piece = NULL;

  if (info->next==NULL)
    return;

  else if (info->next != tail)
    next = info->next->location;

  else
    next = strchr (p, '\0');

  while (p != next)
    {

#pragma omp atomic
      read_num++;

      temp_start = p;

      if (p == '\0' || p == NULL)
	break;

      p = strchr (p, '\n') + 1;

      temp_end = strstr (p, "\n@");

      if (!temp_end)
	temp_end = strchr (p, '\0');

      int result = fastq_read_check (p, strchr (p, '\n') - p, "normal", bl);

      if (result == 1)
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
      else if (result == 0)
	{
#pragma omp critical
	  {
	    //printf("in\n");
	    memcpy (contam, temp_start, temp_end - temp_start);
	    //printf("??\n");
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

  int length = distance;

  int signal = 0, pre_kmer = -1;

  int count_mis = 0, label_m = 0, label_mis = 0, count = 0, match_s = 0;

  char *previous, *key = (char *) malloc (k_mer * sizeof (char) + 1);

//printf("k_mer->%d\n",k_mer);

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
  if (((float) match_s / (float) length) >= tole_rate)
    return 0;
  else
    return 1;
}

/*-------------------------------------*/
void
fasta_process (bloom * bl, Queue * info)
{
  printf ("fasta processing...\n");

  int read_num = 0, read_contam = 0;

  char *p = info->location;

  char *next;

  char *temp = p;

  if (info->next==NULL)
    return;
  else if (info->next != tail)
    next = info->next->location;
  else
    next = strchr (p, '\0');

  while (p != next)
    {
#pragma omp atomic
      read_num++;
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
#pragma omp atomic
	  read_contam++;
#pragma omp critical
	  {
	    memcpy (contam, p, temp - p);
	    contam += (temp - p);
	  }
	}
      p = temp;
    }
  printf ("all->%d\ncontam->%d\n", read_num, read_contam);
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

	  memcpy (temp_key + k_mer - strlen (key), key, sizeof (char) * (strlen (key) + 1));

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
	      return fasta_full_check (bl, begin, next,model);
	    }
	}			//outside if

      else
	{
	  if (!bloom_check (bl, key))
	    {
	      return fasta_full_check (bl, begin, next,model);
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
save_result (char *source, char *obj_file)
{
  printf ("saving...\n");
  char *match = (char *) malloc (400 * sizeof (char)),
    *mismatch = (char *) malloc (400 * sizeof (char)),
    *so_name = (char *) malloc (200 * sizeof (char)),
    *obj_name = (char *) malloc (200 * sizeof (char));

  memset (match, 0, 200);
  memset (mismatch, 0, 200);
  memset (so_name, 0, 200);
  memset (obj_name, 0, 200);

  char *so;
  ((so = strrchr (source, '/'))) && (so += 1, 1) || (so = NULL);

  char *obj;
  ((obj = strrchr (obj_file, '/'))) && (obj += 1, 1) || (obj = NULL);

  if (so)
    strncat (so_name, so, strrchr (source, '.') - so);
  else
    strncat (so_name, source, strrchr (source, '.') - source);

  if (obj)
    strncat (obj_name, obj, strrchr (obj_file, '.') - obj);
  else
    strncat (obj_name, obj_file, strrchr (obj_file, '.') - obj_file);

  if (prefix)
    {
      strcat (match, prefix);
      strcat (mismatch, prefix);
    }
  else if (so)
    {
      strncat (match, source, so - source);
      strncat (mismatch, source, so - source);
    }
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  strcat (match, so_name);
  strcat (mismatch, so_name);
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  strcat (match, "_contam_");
  strcat (mismatch, "_clean_");
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  strcat (match, obj_name);
  strcat (mismatch, obj_name);
  //printf ("match->%s\n", match);
  //printf ("mismatch->%s\n", mismatch);
  if (type == 1)
    {
      strcat (match, ".fasta");
      strcat (mismatch, ".fasta");
    }
  else
    {
      strcat (match, ".fastq");
      strcat (mismatch, ".fastq");
    }
  printf ("match->%s\n", match);
  printf ("mis->%s\n", mismatch);

  write_result (match, contam2);

  write_result (mismatch, clean2);

  free (match);

  free (mismatch);

  free (so_name);

  free (obj_name);

  memset (contam2, 0, strlen (contam2));

  memset (clean2, 0, strlen (clean2));

  clean = clean2;

  contam = contam2;

}
