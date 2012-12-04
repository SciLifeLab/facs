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

static int
build_usage(void)
{
    fprintf(stderr, "\nUsage: ./facs build [options]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-r path/to/something.fasta\n");
    fprintf(stderr, "\t-o have to use it, but does not write file! :_/ XXX\n");
    return 1;
}

int
build_main (int argc, char **argv)
{
  if (argc < 2) return build_usage();

  char *position;
  BIGNUM capacity;

/*-------defaults for bloom filter building-------*/ 
  int opt;
  int k_mer = 21;
  float error_rate = 0.0005;
  float tole_rate = 0.8;

  char* prefix = NULL;
  char* list = NULL;
  char* target_path = NULL;
  char* source = NULL;

  while ((opt = getopt (argc, argv, "ekpo:r:lh")) != -1) {
      switch (opt) {
          case 'e':
              (optarg) && ((error_rate = atof (optarg)), 1);
              break;
          case 'k':
              (optarg) && ((k_mer = atoi (optarg)), 1);
              break;
          case 'p':    
              (optarg) && ((prefix = optarg), 1);
              break;
          case 'o':
              (optarg) && ((target_path = optarg), 1); 
              break;
          case 'r':  
              (optarg) && (source = optarg, 1);  
              break;
          case 'l':
              (optarg) && (list = optarg, 1);  
              break;
          case 'h':
              return build_usage();
          case '?':
              printf ("Unknown option: -%c\n", (char) optopt);
              return build_usage();
      } 
  } 

  build(source, target_path, k_mer, error_rate);

/*
 * @tzcoolman: Is this just replicated code that could be moved to build() function?
 *
  bloom *bl_2 = NEW (bloom);
  Queue *head = NEW (Queue);
  Queue *tail = NEW (Queue);
  head->next = tail;
  F_set *File_head = NEW (F_set);
  File_head = make_list (source, list);
  
  while (File_head) {
      printf ("File_head->%s\n", File_head->filename);
      //map query- into memory--------------
      position = mmaping (File_head->filename);
      if (*position == '>')
    	capacity = strlen (position);
      else
    	capacity = strlen (position) / 2;
      
      init_bloom (bl_2, capacity, error_rate, k_mer);
      ref_add (bl_2, position);
      save_bloom (File_head->filename, bl_2, prefix, target_path);
      bloom_destroy (bl_2);
      
      munmap (position, strlen (position));
      File_head = File_head->next;
    }
*/
  return 0;
}

int
build(char *ref_name, char *target_path, int k_mer, double error_rate)
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

      while (position[bl->k_mer - 1] != '\n') {
    	  memcpy (key, position, sizeof (char) * bl->k_mer);
	      key[bl->k_mer] = '\0';
	      bloom_add (bl, key);
	      position++;
      }

      position += bl->k_mer;
      position = strchr (position, '\n') + 1;
      char *v = strchr (position, '\n');

      if (!v) break;
      else position = v + 1;

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
