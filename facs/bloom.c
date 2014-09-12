#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include "bloom.h"
#include "hashes.h"
#include "file_dir.h"
/*---------------------------*/
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>
/*---------------------------*/
/*
bloom_init, bloom_check, bloom_test, set and set functions are originally from bloom::faster
*/ 
/*---------------------------*/
int seed[20] = {152501029, 152501717,152503097,152500171,152500157,152504837,10161313,10371313,10431313,10501313,10581313,10611313,10641313,10651313,10671313,10731313,10821313,10881313,10951313,11001313};
/*---------------------------*/
int bloom_init (bloom * bloom, BIGNUM size, BIGNUM capacity, double error_rate,
	    int hashes, hash_t hash, int flags)
{
  if (size < 1)
  {
    	fprintf (stderr, "overflow1\n");
    	return -1;
  }
  else
  {
       /* this may waste a little time, but we need to ensure
        * that our array has a prime number of elements, even
        * if we have been requested to do otherwise */
  	bloom->stat.elements = find_close_prime (size);
  }
  if (hashes < 1)
  {
#ifdef DEBUG
  	fprintf (stderr, "hashes was %d,size %lld\n", hashes, size);
#endif
  	return -1;
  }
  else
  {
  	bloom->stat.ideal_hashes = hashes;
  }
  if (hash == NULL)
  {
  	bloom->hash = (hash_t) HASH_FUNC;
  }
  else
  {
  	bloom->hash = hash;
  }
  bloom->inserts = 0;
	/**
	If error rate and capacity were not specified, but size and num hashes were,
	we can calculate the missing elements.
	**/
  if (capacity == 0 || error_rate == 0)
  {
      // From wikipedia, num hashes k that minimizes probability of error is k =~ (0.7 m) / n
      // Therefore n =~ (0.7 m) / k
  	bloom->stat.capacity = 0.7 * bloom->stat.elements / hashes;
  	bloom->stat.e = powf (2.0, (float) -1 * hashes);
  }
  else
  {
  	bloom->stat.capacity = capacity;
  	bloom->stat.e = error_rate;
  }
#ifdef DEBUG
  fprintf (stderr, "bloom_init(%lld,%d) => (%lld,%d) =>%f\n", (BIGCAST) size, hashes, (BIGCAST) bloom->stat.elements, bloom->stat.ideal_hashes, bloom->stat.e);
#endif
  if ((size > TOPLIMIT))
  {
  	fprintf (stderr, "overflow2\n");
  	return -2;
  }
  /* allocate our array of bytes.  where m is the size of our desired 
   * bit vector, we allocate m/8 + 1 bytes. */
  if ((bloom->vector = (char *) calloc (1, sizeof (char) *((BIGNUM) (bloom->stat.elements /8) + 1))) == NULL)
  {
  	perror ("malloc");
  	return -1;
  }
  else
  	memset (bloom->vector, 0, sizeof (char) * ((BIGNUM) (bloom->stat.elements / 8) + 1));
  return 0;
}

void bloom_destroy (bloom * bloom)
{
	free (bloom->vector);
	bloom->vector = NULL;
}

int bloom_check (bloom * bloom, char *str)
{
	int result = 0;
	result = bloom_test (bloom, str, RO);
	return result;
}

void normal_lower(char *str, int length)
{
	char *pstr = str;
	while(length>=0)
	{
		pstr[length] = (char)tolower(pstr[length]);
		length--;
	}
}

int bloom_add (bloom * bloom, char *str)
{
  int ret;
  char* pstr = str;
  //normalize sequence to lowercase
  do
      *pstr = (char)tolower(*pstr);
  while (*pstr++);
  ret = bloom_test (bloom, str, SET);
  if (ret == 0)
  {
      bloom->inserts++;
  }
  return ret;
}

int bloom_test (bloom * bloom, char *str, int mode)
{
  int i, hit=1;
  BIGNUM ret;
  /* as many times as our ideal hash count dictates, salt our key
   * and hash it into the bit vector */
  for (i = 0; i < bloom->stat.ideal_hashes; i++)
  {
      ret = bloom_hash (bloom, str, i, bloom->k_mer);
      if (!test (bloom->vector, ret))
      {
	  hit = 0;
	  if (mode == SET)
	  {
	      set (bloom->vector, ret);
	  }
	  else
	  {
	      /* if we are merely testing, we are done. */
	      return hit;
	  }
      }
  }
  return hit;
}

BIGNUM bloom_hash (bloom * bloom, char *str, int i, int length)
{
	BIGNUM ret = 0;
	ret = (BIGNUM) hash5 (str, seed[i], length) % (BIGNUM) bloom->stat.elements;
	return ret;
}

int set (char *big, BIGNUM index)
{
	deref dr;
	finder (index, &dr);
	big[dr.index] += dr.spot;
	return 0;
}

int test (char *big, BIGNUM index)
{
	deref dr;
	char bucket;
	finder (index, &dr);
	bucket = big[dr.index];
	if ((bucket & dr.spot) == dr.spot)
	{
      		return 1;
  	}
  	else
  	{
      		return 0;
  	}
}

int finder (BIGNUM index, deref * dr)
{

	dr->index = (BIGNUM) (index >> 3);
	dr->spot = 1 << (index & 0x07);
	return 0;
}

char *prefix_make (char *filename, char *prefix, char *target)
{
  char *bloom_file = (char *) calloc (1, 3*ONE * sizeof (char));
  if (target)
  {
  	strcat (bloom_file, target);
  }
  if (!target||is_dir(target))
  {	
	strncat(bloom_file,filename,strrchr(filename,'.')-filename);
	strcat (bloom_file,".bloom");
  }
  return bloom_file;
}

int save_bloom (char *filename, bloom * bl, char *prefix, char *target)
{
  char *bloom_file = NULL;
  int fd;
  bloom_file = prefix_make (filename, prefix, target);
#ifdef DEBUG
  printf ("Bloom file to be written in: %s\n", bloom_file);
#endif
#ifdef __APPLE__
  fd = open (bloom_file, O_RDWR | O_CREAT, PERMS);
#else // assume linux
  #ifndef __clang__
  fd = open (bloom_file, O_RDWR | O_CREAT | O_LARGEFILE, PERMS);
  #endif
#endif
  if (fd < 0)
  {
      perror (bloom_file);
      return -1;
  }
  BIGNUM total_size = sizeof (bloom) + sizeof (char) * ((BIGNUM) (bl->stat.elements / 8) + 1) + sizeof (int) * (bl->stat.ideal_hashes + 1);

#ifdef __APPLE__
  if (ftruncate (fd, total_size) < 0)
#else
  #ifndef __clang__
  if (ftruncate64 (fd, total_size) < 0)
  #endif
#endif
  {
      printf ("[%d]-ftruncate64 error: %s/n", errno, strerror (errno));
      close (fd);
      return 0;
  }
  if (write (fd, bl, sizeof (bloom)) < 0)
  {
      perror (" error writing bloom file ");
      close (fd);
      exit (EXIT_FAILURE);
  }
  total_size = (BIGNUM) (bl->stat.elements / 8) + 1;
  BIGNUM off = 0;
  while (total_size > ONEG*2)
  {
      if (write (fd, bl->vector + off, sizeof (char) * ONEG*2) < 0)
      {
	  perror (" error writing bloom file ");
	  close (fd);
	  exit (EXIT_FAILURE);
      }
      total_size -= ONEG*2;
      off += ONEG*2;
  }
  if (write (fd, bl->vector + off, sizeof (char) * total_size) < 0)
  {
      perror (" error writing bloom file ");
      close (fd);
      exit (EXIT_FAILURE);
  };
  close (fd);

//  memset (bl->vector, 0, sizeof (char) * ((BIGNUM) (bl->stat.elements / 8) + 1));

#ifdef DEBUG
  printf ("big file process OK\n");
#endif
  return 0;

}

int
load_bloom (char *filename, bloom * bl)
{
  int fd;
  int ret;
  BIGNUM off = 0, total_size = 0;
#ifdef DEBUG
  printf ("bloom name->%s\n", filename);
#endif
#ifdef __APPLE__
  fd = open (filename, O_RDONLY, PERMS);
#else
  fd = open64 (filename, O_RDONLY, PERMS);
#endif
  if (fd < 0)
  {
      perror (filename);
      return -1;
  }
  if (read (fd, bl, sizeof (bloom)) < 0)
  {
      perror ("Problem reading Bloom filter");
      close (fd);
      return 0;
  };
  bl->vector = (char *) calloc (1, sizeof (char) * ((BIGNUM) (bl->stat.elements / 8) + 1));
  total_size = ((BIGNUM) (bl->stat.elements / 8) + 1);
  while (total_size > ONEG*2)
  {
  	ret = read (fd, bl->vector + off, sizeof (char) * ONEG*2);
  	if (ret < 0)
	{
		perror ("Problem reading Bloom filter");
		close (fd);
		return 0;
	}
   	total_size -= ONEG*2;
  	off += ONEG*2;
  }
  ret = read (fd, bl->vector + off, sizeof (char) * total_size);
#ifdef DEBUG
  if (ret > 0)
  	printf ("Bloom filter read successfully\n");
  else
  	ret = errno;
#endif
  close (fd);
  return ret;
}

void write_default (char *clean, char *contam, BIGCAST sign)
{
	if (sign==-1)
	{
		fprintf(stdout,"%s",clean);
		fprintf(stderr,"%s",contam);
	}
	else
	{
		fprintf(stdout,"%s",clean);
		fprintf(stderr,"%s",contam);
	}
}
void write_result (char *filename, char *detail)
{
	FILE *fd;
	fd = fopen(filename,"a+");
	if (fd == NULL)
	{
		perror (" error writing result file ");
      		exit (EXIT_FAILURE);
	}
	if (feof(fd)==0)
		fprintf(fd, "%s", detail);  //first time writing data into files normal case
	else
		fprintf(fd, "\n%s",detail); //not first time writing stuff into files, add an extra '\n'
	fclose(fd);
}

void rev_trans (char *s, int length)
{
  int i;
  int j;
  for (i = 0, j = length - 1; i < j; ++i, --j)
  {
      char c = s[i];
      s[i] = s[j];
      s[j] = c;
  }
  i = 0;
  while (i < length)
  {
      switch (s[i])
      {
	case 'A':
	  s[i] = 'T';
	  break;
	case 'C':
	  s[i] = 'G';
	  break;
	case 'G':
	  s[i] = 'C';
	  break;
	case 'T':
	  s[i] = 'A';
	  break;
	case 'a':
	  s[i] = 't';
	  break;
	case 'c':
	  s[i] = 'g';
	  break;
	case 'g':
	  s[i] = 'c';
	  break;
	case 't':
	  s[i] = 'a';
	  break;
	default:
	  break;
      }
      i++;
   }
}

char *mmaping (char *source)
{
  struct stat statbuf;

  int fd = 0;
  char *sm = NULL;

#ifndef __clang__
  if ((fd = open (source, O_RDONLY | O_LARGEFILE)) < 0) {
#else
  if ((fd = open (source, O_RDONLY)) < 0) {
#endif
      fprintf (stderr, "%s: %s\n", source, strerror (errno));
      exit (EXIT_FAILURE);
  }
  if (fstat (fd, &statbuf) < 0) {
      fprintf (stderr, "%s: %s\n", source, strerror (errno));
      close (fd);
      exit (EXIT_FAILURE);
  } else if (statbuf.st_size == 0) {
      fprintf (stderr, "%s: %s\n", source, "File is empty");
      close (fd);
      exit (EXIT_FAILURE);
  }

  sm = mmap (0, (BIGCAST) statbuf.st_size, PROT_READ, MAP_SHARED | MAP_NORESERVE, fd, 0);

  if (MAP_FAILED == sm){
      fprintf (stderr, "%s: %s\n", source, strerror (errno));
      close (fd);
      exit (EXIT_FAILURE);
  }

  close (fd);
  return sm;
}

char *large_load (char *fifoname)
{
  int x = 0;
  char ch;
  FILE *fd;
#ifdef __APPLE__
  fd = fopen (fifoname, "r");
#else
  fd = fopen64 (fifoname, "r");
#endif
  char *data = (char *) calloc (1, (ONEG+1) * sizeof (char));
  data[ONEG] = '\0';
  while ((ch = fgetc (fd)) != EOF)
  {
      data[x] = ch;
      x++;
  }
#ifdef DEBUG
  printf ("data length->%lld\n", (BIGNUM) strlen (data));
#endif
  fclose (fd);
  return data;
}

BIGCAST get_size (char *filename)
{
  struct stat buf;
  if (stat (filename, &buf) != -1)
  {
  	return buf.st_size;
  }
  else
  {
  	return 0;
  }
}

void print_bloom_info(bloom *bl) {
   printf("Word size:\t%d\n", bl->k_mer);
   printf("Filter size:\t%lld\n", bl->stat.elements);
   printf("Max error rate:\t%e\n", bl->stat.e);
   printf("#inserts:\t%lld\n", bl->inserts);
 }
