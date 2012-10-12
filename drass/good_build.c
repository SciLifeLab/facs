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

#define PERMS 0600
#define NEW(type) (type *) malloc(sizeof(type))
#define FILE_MODE (S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)
/*-------------------------------------*/
long sec, usec, i;
struct timezone tz;
struct timeval tv, tv2;
struct stat statbuf;
/*-------------------------------------*/
float error_rate, tole_rate, contamination_rate;
BIGNUM total_size = 0;
int k_mer = 21, mode = 0, count = 0, type;
long long hit = 0, un_hit = 0;

BIGNUM capacity;

char *source, *position, *prefix, *list;

Queue *head, *tail, *head2;

F_set *File_head, *File_head2;

bloom *bl_2;

void struc_init ();
void init_bloom (bloom * bl);
void init_query (char *source);
void init (int argc, char **argv);
void ref_add (bloom * bl, char *position);
void fastq_add (bloom * bl, char *position);
void fasta_add (bloom * bl, char *position);

char *fasta_data (bloom * bl_2, char *data);
int build(char *ref_name, char *target_path, int k_mer, double error_rate);


int main (int argc, char *argv[])
{
  //gettimeofday (&tv, &tz);
  build("k_12.fasta", argv[0], 21, 0.0005);
/*
  init (argc, argv);
  struc_init ();

  while (File_head)
    {
      printf ("File_head->%s\n", File_head->filename);

      init_query (File_head->filename);
      init_bloom (bl_2);

      ref_add (bl_2, position);
      save_bloom (File_head->filename, bl_2, prefix, argv[0]);

      bloom_destroy (bl_2);

      if (!strstr (File_head->filename, ".fifo"))
	munmap (position, strlen (position));
      else
	free (position);
      File_head = File_head->next;
    }

  gettimeofday (&tv2, &tz);

  sec = tv2.tv_sec - tv.tv_sec;
  usec = tv2.tv_usec - tv.tv_usec;

#ifdef DEBUG
  printf ("all finished...\n");
  printf ("total=%ld sec\n", sec/100000000);
  printf ("Same K_mer->%ld\n,New K_mer->%ld\n", hit, un_hit);
#endif
*/
  return 0;
}

/*-------------------------------------*/
int build(char *ref_name, char *target_path, int k_mer, double error_rate)
{
    printf("ERROR RATE, just coming from python iface: %d %f\n", k_mer, error_rate);
	char *position = mmaping (ref_name);
	
	bloom *bl = NEW (bloom);
	bl->k_mer = k_mer;
    bl->stat.e = error_rate;

    bl->stat.capacity = strlen(position);
    get_rec(&bl->stat);

	bloom_init(bl, bl->stat.elements, bl->stat.capacity, bl->stat.e,bl->stat.ideal_hashes, NULL, 3);
	ref_add (bl, position);
	save_bloom (ref_name, bl, NULL, target_path);
	
	return 0;
}

void
struc_init ()
{
  bl_2 = NEW (bloom);
  head = NEW (Queue);
  tail = NEW (Queue);
  head->next = tail;
  head2 = head;
  File_head = NEW (F_set);
  File_head = make_list (source, list);
  File_head = File_head->next;
}

/*-------------------------------------*/
void
init_query (char *source)
{
  if (strstr (source, ".fifo"))
    position = large_load (source);
  else
    position = mmaping (source);

  if (*position == '>')
    capacity = strlen (position);
  else
    capacity = strlen (position) / 2;
}

/*-------------------------------------*/

void
init (int argc, char **argv)
{
  if (argc == 1 || !strcmp (argv[1], "-h") || !strcmp (argv[1], "-help"))
    {
      help ();
      build_help ();
    }
/*-------default-------*/
  type = 2;
  mode = 1;
  k_mer = 21;
  error_rate = 0.0005;
  prefix = NULL;
/*-------default-------*/
  int x;
  while ((x = getopt (argc, argv, "e:k:m:o:r:l:h")) != -1)
    {
      switch (x)
	{
	case 'e':
	  (optarg) && ((error_rate = strtod (optarg, NULL)));
	  break;
	case 'k':
	  (optarg) && ((k_mer = atoi (optarg)), 1);
	  break;
	case 'm':
	  (optarg) && ((mode = atoi (optarg)), 1);
	  if (mode > 2 || mode < 1)
	    {
	      perror ("Mode select error.");
	      exit (0);
	    }
	  break;
	case 'o':
	  (optarg) && ((prefix = optarg), 1);
	  break;
	case 'r':
	  (optarg) && ((source = optarg), 1);
	  break;
	case 'l':
	  (optarg) && (list = optarg, 1);
	  break;
	case 'h':
	  help ();
	  remove_help ();
	  break;
	case '?':
	  printf ("Unknown option: -%c\n", (char) optopt);
	  exit (0);
	}

    }

  if ((!source) && (!list))
    {
      perror ("No source.");
      exit (0);
    }
}

/*-------------------------------------*/
void
init_bloom (bloom * bl)
{
  int flags = 3;
  hash_t hash = NULL;

  get_suggestion(&bl->stat, capacity, error_rate);

#ifdef DEBUG
  printf ("Capacity: %lld\n", bl->stat.capacity);
  printf ("Vector size: %lld\n", bl->stat.elements);
  printf ("Ideal hashes: %d\n", bl->stat.ideal_hashes);
  printf ("Error rate: %f\n", bl->stat.e);
  printf ("Real size: %lld\n", bl->stat.elements / 8);
#endif

  bloom_init (bl, bl->stat.elements, bl->stat.capacity, bl->stat.e,
	      bl->stat.ideal_hashes, hash, flags);

  bl->k_mer = k_mer;

}

/*-------------------------------------*/
char *
fasta_title (char *full)
{

  char *ptr;

  ptr = strchr (full, '\n');

  return ptr + 1;
}

/*-------------------------------------*/
void
fasta_add (bloom * bl, char *position)
{
  while (*position != '\0')
    {
      if (*position == '>')
	position = fasta_title(position);
      else
	position = fasta_data(bl, position);
    }
}

/*-------------------------------------*/
void
fastq_add (bloom * bl, char *position)
{

  char *key = (char *) malloc (bl->k_mer * sizeof (char) + 1);

  while (position[0] != '\0')
    {
      position = strchr (position, '\n') + 1;

      while (position[bl->k_mer - 1] != '\n')
	{
	  memcpy (key, position, sizeof (char) * bl->k_mer);

	  key[bl->k_mer] = '\0';

	  if (bloom_add (bl, key))
	    hit++;
	  else
	    un_hit++;

	  position++;
	}

      position += bl->k_mer;

      position = strchr (position, '\n') + 1;

      char *v = strchr (position, '\n');

      if (!v)
	break;
      else

	position = v + 1;

    }
  free (key);
}

/*-------------------------------------*/
char *
fasta_data (bloom * bl_2, char *data)
{

  char *key = (char *) malloc (bl_2->k_mer * sizeof (char) + 1);
  char *p = data;
  int n = 0, m = 0;

  while (*p != '>' && *p != '\0') {
      while (n < bl_2->k_mer) {
	    if (p[m] == '>' || p[m] == '\0') {
	      m--;
	      break;
	    }

	  if (p[m] != '\r' && p[m] != '\n')
	    key[n++] = p[m];
	    m++;
      }

      key[n] = '\0';


      if (strlen(key) == bl_2->k_mer) {

	  if (bloom_add (bl_2, key))
	    hit++;
	  else
	    un_hit++;
      }

      p += 1;
      n = 0;
      m = 0;
      }

  free(key);
  return p;
}

/*-------------------------------------*/
void
ref_add (bloom * bl, char *position)
{
  if (*position == '>')
    fasta_add (bl, position);
  else if (*position == '@')
    fastq_add (bl, position);
  else
    {
      perror ("wrong format\n");
      exit (-1);
    }
}
