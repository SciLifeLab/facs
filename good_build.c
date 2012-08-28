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
int k_mer = 21, mode = 0, count = 0;
long long hit = 0, un_hit = 0;

BIGNUM capacity;

char *dm, *dst, *mis_dm, *mis_dst, *source,
     *obj, *temp_title, *position, *prefix;

Queue *head, *tail;

void init_bloom (bloom * bl);
void init (int argc, char **argv);
void fastq_add (bloom * bl, char *position);
void fasta_add (bloom * bl, char *position);

char *mmaping (char *filename);
char *fasta_title (char *full);
char *large_load (char *filename);
char *fasta_data (bloom * bl_2, char *data);

long long get_size (char *strFileName);

main (int argc, char *argv[])
{
  gettimeofday (&tv, &tz);	// time test
  bloom *bl_2;
  bl_2 = NEW (bloom);
  head = NEW (Queue);
  tail = NEW (Queue);
  head->next = tail;
 
  init (argc, argv);

  char *program_path = (char *) malloc (200 * sizeof (char));

  position = mmaping(source);
  
  if (*position == '>')
    capacity = strlen (position);
  else
    capacity = strlen (position) / 2;
  
  init_bloom (bl_2);

  init_bloom (bl_2);
  fasta_add(bl_2, position);

#ifdef DEBUG
  strcat (program_path, argv[0]);
  printf ("program_path->%s\n",program_path);
#endif

  save_bloom (source, bl_2, prefix, argv[0]);
  munmap (position, statbuf.st_size);

  gettimeofday (&tv2, &tz);

  sec = tv2.tv_sec - tv.tv_sec;
  usec = tv2.tv_usec - tv.tv_usec;

#ifdef DEBUG
  printf ("all finished...\n");
  printf ("total=%ld sec\n", sec);
  printf ("Same K_mer->%ld\n,New K_mer->%ld\n", hit, un_hit);
#endif

  return 0;
}

/*-------------------------------------*/
char *
mmaping (char *source)
{

  int src;
  char *sm;


#ifdef __APPLE__
  src = open (source, O_RDONLY,0644);
#else
  src = open (source, O_RDONLY|O_LARGEFILE,0644);
#endif
  if (src < 0) {
      perror (" open source ");
      exit (EXIT_FAILURE);
  }

  if (fstat (src, &statbuf) < 0) {
      perror (" fstat source ");
      exit (EXIT_FAILURE);
  }

  sm = mmap (0, (long long) statbuf.st_size, PROT_READ,
	   MAP_SHARED | MAP_NORESERVE, src, 0);

  if (MAP_FAILED == sm) {
      perror (" mmap source ");
      exit (EXIT_FAILURE);
  }

  return sm;
}

/*-------------------------------------*/
char *
large_load (char *filename)
{

  int fd;

#ifdef DEBUG
  printf ("queryname->%s\n", filename);
#endif

#ifdef __APPLE__
  fd = open(filename, O_RDWR|O_CREAT, 0644); 
#else
  fd = open(filename, O_RDWR|O_CREAT|O_LARGEFILE, 0644); 
#endif

  if (fd < 0) {
      perror (filename);
      return NULL;
  }

  (strstr (filename, ".fasta") || strstr (filename, ".fna"))
    && ((total_size = get_size (filename)), 1)
    || (total_size = (get_size (filename) * 2));


#ifdef __APPLE__
  if (ftruncate(fd, total_size) < 0)
#else
  if (ftruncate64(fd, total_size) < 0)
#endif
    {
      printf ("[%d]-ftruncate64 error: %s/n", errno, strerror (errno));
      close (fd);
      return 0;
    }

  printf ("total_size->%lld\n", total_size);

  char *data = (char *) malloc ((total_size + 1) * sizeof (char));

  read (fd, data, total_size);
  close (fd);

#ifdef DEBUG
// Too verbose
// printf("data->%s\n",data);
#endif

  return data;
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
  mode = 1;
  k_mer = 21;
  error_rate = 0.0005;
  prefix = NULL;
/*-------default-------*/
  int x;
  while ((x = getopt (argc, argv, "e:k:m:o:r:")) != -1)
    {
      switch (x)
	{
	case 'e':
	  (optarg) && ((error_rate = strtod(optarg, NULL)));
	  break;
	case 'k':
	  (optarg) && ((k_mer = atoi (optarg)), 1);
	  break;
	case 'm':
	  (optarg) && ((mode = atoi (optarg)), 1);
	  if (mode != 1 && mode != 2)
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
	  if (!source)
	    {
	      perror ("error, no resource.");
	      exit (0);
	    }
	  break;
	case '?':
	  printf("Unknown option: -%c\n",(char)optopt);
	  exit (0);
	}

  }

}

/*-------------------------------------*/
void
init_bloom (bloom * bl)
{
  BIGNUM size;

  int status, hashes, flags;

  hash_t hash = NULL;

  flags = 3;

  get_suggestion(&bl->stat, capacity, error_rate);

#ifdef DEBUG
  printf ("Capacity: %lld\n", bl->stat.capacity);
  printf ("Vector size: %lld\n", bl->stat.elements);
  printf ("Ideal hashes: %d\n", bl->stat.ideal_hashes);
  printf ("Error rate: %f\n", bl->stat.e);
  printf ("Real size: %lld\n", bl->stat.elements / 8);
#endif

  bloom_init(bl, bl->stat.elements, bl->stat.capacity, bl->stat.e,
	         bl->stat.ideal_hashes, hash, flags);

  bl->k_mer = k_mer;

}

/*-------------------------------------*/
char *
fasta_title (char *full)
{

  char *ptr;

  ptr = strchr (full, '\n');

  temp_title = (char *) malloc ((ptr - full + 1) * sizeof (char));

  strncpy (temp_title, full, ptr - full + 1);	//include '\n'

  return ptr + 1;
}

/*-------------------------------------*/
void
fasta_add (bloom * bl, char *position)
{
  while (strlen (position) > k_mer + 2)
    {
      if (position[0] == '>')
	position = fasta_title (position);
      else
	position = fasta_data (bl, position);
    }
}

/*-------------------------------------*/
void
fastq_add (bloom * bl, char *position)
{

  char *key = (char *) malloc (k_mer * sizeof (char) + 1);
  char *position2;

  while (position[0] != '\0')
    {
      position = strchr (position, '\n') + 1;

      while (position[k_mer - 1] != '\n')
	{
	  memcpy (key, position, sizeof (char) * k_mer);

	  key[k_mer] = '\0';

	  if (bloom_add (bl, key))
	    hit++;
	  else
	    un_hit++;

	  position++;
	}

      position += k_mer;

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
long long
get_size (char *strFileName)
{
  struct stat temp;
  stat (strFileName, &temp);
  if (strstr (strFileName, ".fna") || strstr (strFileName, ".fasta"))
    return (temp.st_size);
  else
    return (temp.st_size / 2);

}

/*-------------------------------------*/
char *
fasta_data (bloom * bl_2, char *data)
{

  char *key = (char *) malloc (k_mer * sizeof (char) + 1);

  char *p = data;

  int n = 0, m = 0;

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

      if (strlen (key) == k_mer)
	{

	  if (bloom_add (bl_2, key))
	    hit++;
	  else
	    un_hit++;
	}
      p += 1;
      n = 0;
      m = 0;
    }
  free (key);
  return p;
}

/*-------------------------------------*/
