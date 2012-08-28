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
/*-------------------------------------*/
float error_rate, sampling_rate, contamination_rate, tole_rate;
/*-------------------------------------*/
int k_mer = 0, mode, mytask, ntask, type = 2, reads_num =
  0, reads_contam = 0, checky = 0;
/*-------------------------------------*/
char *source, *all_ref, *position, *prefix, *clean, *contam, *clean2,
  *contam2, *query_list;
/*-------------------------------------*/
Queue *head, *tail, *head2, *tail2;
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
void evaluate (char *detail, char *filename);
void save_result (char *source, char *obj_file);
void statistic_save (char *detail, char *filename);
/*-------------------------------------*/
int fastq_full_check (bloom * bl, char *p, int distance);
int fasta_full_check (bloom * bl, char *begin, char *next, char *model);
int fastq_read_check (char *begin, int length, char *model, bloom * bl);
int fasta_read_check (char *begin, char *next, char *model, bloom * bl);
/*-------------------------------------*/
char *jump (char *target);
char *mmaping (char *source);
/*-------------------------------------*/
//long long get_size(char * strFileName);
/*-------------------------------------*/
//bloom *load_bloom (char *filename, bloom * bl);
/*-------------------------------------*/
main (int argc, char **argv)
{

  long sec, usec, i;
  struct timezone tz;
  struct timeval tv, tv2;

  gettimeofday (&tv, &tz);

  init (argc, argv);		//initialize 
  struc_init ();

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
    }				// End of single - no implied barrier (nowait)
  }				// End of parallel region - implied barrier

#ifdef DEBUG
  printf ("finish processing...\n");
#endif

  evaluate (detail, source);
  statistic_save (detail, source);
  munmap (position, statbuf.st_size);
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
      check_help ();
    }
/*-------default-------*/
  mode = 1;
  k_mer = 21;
  tole_rate = 0.8;
  error_rate = 0.0005;
  sampling_rate = 1;
  prefix = NULL;
/*-------default-------*/

  int x;
  while ((x = getopt (argc, argv, "e:k:m:t:o:r:q:s:")) != -1)
    {
      //printf("optind: %d\n", optind);
      switch (x)
	{
	case 'e':
	  (optarg) && ((error_rate = atof (optarg)), 1);
	  break;
	case 'k':
	  (optarg) && ((k_mer = atoi (optarg)), 1);
	  break;
	case 'm':
	  (optarg) && ((mode = atoi (optarg)), 1);
	  break;
	case 't':
	  (optarg) && ((tole_rate = atof (optarg)), 1);
	  break;
	case 's':
	  (optarg) && ((sampling_rate = atof (optarg)), 1);
	  break;
	case 'o':
	  (optarg) && ((prefix = optarg), 1);
	  break;
	case 'r':
	  (optarg) && ((all_ref = optarg), 1);
	  break;
	case 'q':
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
    {
      perror ("No source.");
      exit (0);
    }

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
  head = NEW (Queue);
  tail = NEW (Queue);
  head->next = tail;
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

  if (type == 1) {
      for (add = 0; add < cores; add++)
	{
	  Queue *x = NEW (Queue);

          if (add == 0 && *full != '>')

	  temp = strchr (full, '>');	//drop the possible fragment

	  if (add != 0)
	    temp = strchr (full + offsett, '>');


	  x->location = temp;
	  x->number = add;
	  x->next = pos->next;
	  pos->next = x;
	  pos = pos->next;
	}
  } else {

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
	} 

   }

  return;
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

  return sm;
}

/*-------------------------------------*/
void
fastq_process (bloom * bl, Queue * info)
{

#ifdef DEBUG
  printf("fastq processing...\n");
#endif

  char *p = info->location;
  char *next, *temp, *temp_piece = NULL;

  if (info->next==NULL)
    return;

  else if (info->next != tail)
    next = info->next->location;

  else
    next = strchr (p, '\0');

  while (p != next)
    {
      temp = jump (p);		//generate random number and judge if need to scan this read

      if (p != temp)
	{
	  p = temp;
	  continue;
	}

      if (*p == '\0' || !p)
	break;

#pragma omp atomic
      reads_num++;

      p = strchr (p, '\n') + 1;
      int distance = strchr (p, '\n') - p;

      if (!fastq_read_check (p, distance, "normal", bl))
#pragma omp atomic
	reads_contam++;

      p = strchr (p, '\n') + 1;
      p = strchr (p, '\n') + 1;
      p = strchr (p, '\n') + 1;
    }				// outside while
  if (temp_piece)
    free (temp_piece);
}

/*-------------------------------------*/
int
fastq_read_check (char *begin, int length, char *model, bloom * bl)
{
  char *p = begin;
  int distance = length, label_m = 0, label_mis = 0;
  int signal = 0, pre_kmer = 10, count_mis = 0;
  char *previous, *key = (char *) malloc (k_mer * sizeof (char) + 1);

  while (distance > 0)
    {
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

          if (mode == 1) {
              if (bloom_check (bl, key))
                  return fastq_full_check (bl, begin, length);
	  } else {

          if (!bloom_check (bl, key)) {
              return fastq_full_check (bl, begin, length);
          }
	}
  }				// inner while

  free (key);

  if (model == "normal")	//use recursion to check the sequence forward and backward
    return fastq_read_check (begin, length, "reverse", bl);
  else
    {
      if (mode == 1)
	label_mis += distance;
      else
	label_m += distance;
      return 1;
    }
}

/*-------------------------------------*/
int
fastq_full_check (bloom * bl, char *p, int distance)
{

  //printf("fastq full check...\n");

  int signal = 0, pre_kmer = -1;

  int count_mis = 0, label_m = 0, label_mis = 0, count = 0, match_s = 0;

  char *previous, *key = (char *) malloc (k_mer * sizeof (char) + 1);

#pragma omp atomic
  checky++;

  int length = distance;

  while (distance >= k_mer)
    {
      memcpy (key, p, sizeof (char) * k_mer);
      key[k_mer] = '\0';
      previous = p;
      p += 1;

      if (bloom_check (bl, key)) {

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
        } else {
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

#ifdef DEBUG
  printf("fasta processing...\n");
#endif

  char *p = info->location;
  char *temp_next, *next, *temp, *temp_piece = NULL;

  if (info->next==NULL)
    return;
  else if (info->next != tail)
    next = info->next->location;
  else
    next = strchr (p, '\0');

  while (p != next)
    {

      temp = jump (p);		//generate random number and judge if need to scan this read

      if (p != temp)
	{
	  p = temp;
	  continue;
	}

#pragma omp atomic
      reads_num++;

      temp_next = strchr (p + 1, '>');
      if (!temp_next)
	temp_next = next;

      if (!fasta_read_check (p, temp_next, "normal", bl))
	{
#pragma omp atomic
	  reads_contam++;
	}

      p = temp_next;
    }
}

/*-------------------------------------*/
int
fasta_read_check (char *begin, char *next, char *model, bloom * bl)
{

//printf("fasta read check...\n");

  char *p = strchr (begin + 1, '\n') + 1;

  if (!p || *p == '>')
     return 1;

  char *start = p;

  char *key = (char *) malloc ((k_mer + 1) * sizeof (char));

  char *pre_key = (char *) malloc ((k_mer + 1) * sizeof (char));

  int n = 0, m = 0, count_enter = 0, result = 0;

  int label_m = 0, label_mis = 0;

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

  //printf("fasta full check...\n");

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
evaluate (char *detail, char *filename)
{
  char buffer[100] = { 0 };

  printf ("all->%d\n", reads_num);

  printf ("contam->%d\n", reads_contam);

//printf("possbile->%d\n",checky);

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

  sprintf (buffer, "all->%d ", reads_num);
  strcat (detail, buffer);

  memset (buffer, 0, 100);
  sprintf (buffer, "contam->%d ", reads_contam);

  strcat (detail, buffer);
  memset (buffer, 0, 100);

  sprintf (buffer, "possbile->%d ", checky);
  strcat (detail, buffer);

  memset (buffer, 0, 100);
  sprintf (buffer, "contamination rate->%f", contamination_rate);

  strcat (detail, buffer);
  memset (buffer, 0, 100);

  reads_num = 0;
  reads_contam = 0;
  checky = 0;
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
