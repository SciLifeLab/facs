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
#include "file_dir.h"
#include "tool.h"
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
  *contam2, *list;
/*-------------------------------------*/
Queue *head, *tail, *head2;
/*-------------------------------------*/
bloom *bl_2;
/*-------------------------------------*/
F_set *File_head;
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

  if (strstr (source, ".fifo"))
    position = large_load (source);
  else
    position = mmaping (source);

  get_parainfo (position);

  char *detail = (char *) malloc (1000 * 1000 * sizeof (char));

  memset (detail, 0, 1000 * 1000);

  while (File_head)
    {

      load_bloom (File_head->filename, bl_2);

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

//if (strstr(source,".fifo"))
//  {
      //source = (char *) malloc (500 * sizeof (char));
//  strncpy(source,fifoname,(strrchr(fifoname,'.')-fifoname));
//  }

      evaluate (detail, File_head->filename);

      File_head = File_head->next;

      head = head2;

      bloom_destroy (bl_2);
    }				//end while

  statistic_save (detail, source);

  if (!strstr (source, ".fifo"))
    munmap (position, strlen (position));

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
  while ((x = getopt (argc, argv, "e:k:m:t:o:r:q:s:l:h")) != -1)
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
	case 'l':
	  (optarg) && (list = optarg, 1);
	  break;
	case 'h':
	  help ();
	  check_help ();
	  break;
	case '?':
	  printf ("Unknown option: -%c\n", (char) optopt);
	  exit (0);
	}
    }

  if (((!all_ref) && (!list)) || (!source))
    {
      perror ("No source.");
      exit (0);
    }

  //if (strstr (source, ".fasta") || strstr (source, ".fna"))
  //type = 1;

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
  head2 = head;
  File_head = NEW (F_set);
  File_head = make_list (all_ref, list);
  File_head = File_head->next;
}

/*-------------------------------------*/
void
get_parainfo (char *full)
{
  printf ("distributing...\n");

  char *temp = full;
  int cores = omp_get_num_procs ();
  int offsett = strlen (full) / cores;
  int add = 0;

  printf ("task->%d\n", offsett);

  Queue *pos = head;

  if (*full == '>')
    type = 1;
  else if (*full == '@')
    type = 2;
  else
    {
      perror ("wrong format\n");
      exit (-1);
    }
  if (type == 1)
    {
      for (add = 0; add < cores; add++)
	{
	  Queue *x = NEW (Queue);

	  if (add == 0 && *full != '>')

	    temp = strchr (full, '>');	//drop the possible fragment

	  if (add != 0)
	    temp = strchr (full + offsett * add, '>');


	  x->location = temp;
	  x->number = add;
	  x->next = pos->next;
	  pos->next = x;
	  pos = pos->next;
	}
    }
  else
    {

      for (add = 0; add < cores; add++)
	{
	  Queue *x = NEW (Queue);

	  if (add == 0 && *full != '@')

	    temp = strstr (full, "\n@") + 1;	//drop the fragment

	  //printf("offset->%d\n",offsett*add);

	  if (add != 0)

	    temp = strstr (full + offsett * add, "\n@");

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
void
fastq_process (bloom * bl, Queue * info)
{

#ifdef DEBUG
  printf ("fastq processing...\n");
#endif

  char *p = info->location;
  char *next, *temp, *temp_piece = NULL;

  if (info->next == NULL)
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

      p = strchr (p, '\n') + 1;
      if (fastq_read_check (p, strchr (p, '\n') - p, "normal", bl, tole_rate)>0)
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
void
fasta_process (bloom * bl, Queue * info)
{

#ifdef DEBUG
  printf ("fasta processing...\n");
#endif

  char *p = info->location;
  char *temp_next, *next, *temp, *temp_piece = NULL;

  if (info->next == NULL)
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

      if (fasta_read_check (p, temp_next, "normal", bl, tole_rate)>0)
	{
#pragma omp atomic
	  reads_contam++;
	}

      p = temp_next;
    }
}
/*-------------------------------------*/
void
evaluate (char *detail, char *filename)
{
  char buffer[200] = { 0 };

  printf ("all->%d\n", reads_num);

  printf ("contam->%d\n", reads_contam);

  printf ("bloomname->%s\n", filename);

  contamination_rate = (double) (reads_contam) / (double) (reads_num);

  if ((mode == 1 && contamination_rate == 0)
      || (mode == 2 && contamination_rate == 1))
    printf ("clean data...\n");

  else if (mode == 1)
    printf ("contamination rate->%f\n", contamination_rate);
  else
    printf ("contamination rate->%f\n", 1 - contamination_rate);

  strcat (detail, "\nxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  strcat (detail, "Bloomfile\tAll\tContam\tcontam_rate\n");

  strcat (detail, filename);
  //strcat (detail, "  ");
  //strcat (detail, "   \n");

  sprintf (buffer, "  %d\t%d\t%f\n", reads_num, reads_contam,
	   contamination_rate);
  strcat (detail, buffer);

  //memset (buffer, 0, 100);
  //sprintf (buffer, "contam->%d ", reads_contam);

  //strcat (detail, buffer);
  //memset (buffer, 0, 100);

  //sprintf (buffer, "possbile->%d ", checky);
  //strcat (detail, buffer);

  //memset (buffer, 0, 100);
  //sprintf (buffer, "contamination rate->%f", contamination_rate);

  //strcat (detail, buffer);
  //memset (buffer, 0, 100);

  reads_num = 0;
  reads_contam = 0;
  checky = 0;
  contamination_rate = 0;

  //printf("detail->%s\n",detail);
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

  printf ("Info name->%s\n", save_file);

  write_result (save_file, detail);
}
