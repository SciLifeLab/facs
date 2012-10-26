#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
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
#include "build.h"
#include "bloom.h"
#include "file_dir.h"
/*
void struc_init ();
void init_bloom (bloom * bl);
void init_query (char *source);
void init (int argc, char **argv);
void ref_add (bloom * bl, char *position);
void fastq_add (bloom * bl, char *position);
void fasta_add (bloom * bl, char *position);

char *fasta_data (bloom * bl_2, char *data);
int build(char *ref_name, char *target_path, int k_mer, double error_rate);
*/

int
build_main (int k_mer, float error_rate, char *source, char *list,
	    char *prefix, char *target_path, int help)
{
  if (help == 1)
    {
      build_help ();
      exit (1);
    }
/*-------------------------------------*/
  char *position;
  BIGNUM capacity;
  //BIGCAST hit = 0, un_hit = 0;
/*-------------------------------------*/
  if (help == 3)
    build ("k_12.fasta", NULL, 21, 0.0005);
/*-------------------------------------*/
  bloom *bl_2 = NEW (bloom);
  Queue *head = NEW (Queue);
  Queue *tail = NEW (Queue);
  head->next = tail;
  F_set *File_head = NEW (F_set);
  File_head = make_list (source, list);
/*-------------------------------------*/
  while (File_head)
    {
      printf ("File_head->%s\n", File_head->filename);
      /*map query- into memory-------------- */
      position = mmaping (File_head->filename);
      if (*position == '>')
	capacity = strlen (position);
      else
	capacity = strlen (position) / 2;
      /*init bloom-------------------------- */
      init_bloom (bl_2, capacity, error_rate, k_mer);
      /*hashing and bloom building---------- */
      ref_add (bl_2, position);
      /*save bloom-------------------------- */
      save_bloom (File_head->filename, bl_2, prefix, target_path);
      /*clear memory------------------------ */
      bloom_destroy (bl_2);
      munmap (position, strlen (position));
      File_head = File_head->next;
    }
  printf ("all finished...\n");
#ifdef DEBUG
  /*int sec, usec;
  struct timeval tv1, tv2;
  struct timezone tz;
  gettimeofday (&tv2, &tz);
  sec = tv2.tv_sec - tv.tv_sec;
  usec = tv2.tv_usec - tv.tv_usec;
  printf ("total=%ld sec\n", sec);*/
  //printf ("Same K_mer->%ld\n,New K_mer->%ld\n", hit, un_hit);
#endif

  return 0;
}

/*-------------------------------------*/
int
build (char *ref_name, char *target_path, int k_mer, double error_rate)
{
  char *position = mmaping (ref_name);

  bloom *bl = NEW (bloom);
  bl->k_mer = k_mer;
  bl->stat.e = error_rate;
  bl->stat.capacity = strlen (position);
  get_rec (&bl->stat);

  bloom_init (bl, bl->stat.elements, bl->stat.capacity, bl->stat.e,
	      bl->stat.ideal_hashes, NULL, 3);
  ref_add (bl, position);
  save_bloom (ref_name, bl, NULL, target_path);

  return 0;
}

/*-------------------------------------*/
void
init_bloom (bloom * bl, BIGNUM capacity, float error_rate, int k_mer)
{
  int flags = 3;

  get_suggestion (&bl->stat, capacity, error_rate);

#ifdef DEBUG
  printf ("Capacity: %lld\n", bl->stat.capacity);
  printf ("Vector size: %lld\n", bl->stat.elements);
  printf ("Ideal hashes: %d\n", bl->stat.ideal_hashes);
  printf ("Error rate: %f\n", bl->stat.e);
  printf ("Real size: %lld\n", bl->stat.elements / 8);
#endif

  bloom_init (bl, bl->stat.elements, bl->stat.capacity, bl->stat.e,
	      bl->stat.ideal_hashes, NULL, flags);

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
	position = fasta_title (position);
      else
	position = fasta_data (bl, position);
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
/*
	  if (bloom_add (bl, key))
	    hit++;
	  else
	    un_hit++;
*/
	  bloom_add (bl, key);

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

  while (*p != '>' && *p != '\0')
    {
      while (n < bl_2->k_mer)
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


      if (strlen (key) == bl_2->k_mer)
	{
/*
	  if (bloom_add (bl_2, key))
	    hit++;
	  else
	    un_hit++;
*/
      }

	  bloom_add (bl_2, key);
	  p += 1;
	  n = 0;
	  m = 0;
	}

      free (key);
      return p;
    }
/*-------------------------------------*/
  void ref_add (bloom * bl, char *position)
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
