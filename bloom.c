#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include "bloom.h"
#include "hashes.h"
/*---------------------------*/
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>
#define TWOG 2000000000
#define hashsize(n) ((BIGNUM)1<<(n))
#define hashmask(n) (hashsize(n) - 1)
#define NEW(type) (type *) malloc(sizeof(type))
#define FILE_MODE (S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)

int seed[20] =
  { 152501029, 152501717, 152503097, 152500171, 152500157, 152504837,
  10161313, 10371313, 10431313, 10501313, 10581313, 10611313, 10641313,
    10651313,
  10671313, 10731313, 10821313, 10881313, 10951313, 11001313
};

int
bloom_init (bloom * bloom, BIGNUM size, BIGNUM capacity, float error_rate,
	    int hashes, hash_t hash, int flags)
{
  if (size < 1)
    {
      fprintf (stderr, "overflow1\n");
      return -1;
    } else {
      /* this may waste a little time, but we need to ensure
       * that our array has a prime number of elements, even
       * if we have been requested to do otherwise */
      bloom->stat.elements = find_close_prime (size);
    }

  if (hashes < 1)
    {
      fprintf (stderr, "hashes was %d,size %lld\n", hashes, size);
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
  if ((flags & BVERBOSE) == BVERBOSE)
    {
      fprintf (stderr, "bloom_init(%lld,%d) => (%lld,%d) =>%f\n",
	       (BIGCAST) size, hashes, (BIGCAST) bloom->stat.elements,
	       bloom->stat.ideal_hashes, bloom->stat.e);
      //verbose = 1;
    }
  else
    {
      //verbose = 0;
    }
  if ((size > TOPLIMIT))
    {
      fprintf (stderr, "overflow2\n");
      return -2;
    }

  /* allocate our array of bytes.  where m is the size of our desired 
   * bit vector, we allocate m/8 + 1 bytes. */
  if ((bloom->vector =
       (char *) malloc (sizeof (char) *
			((long long) (bloom->stat.elements / 8) + 1))) ==
      NULL)
    {
      //printf("-->%d\n",strlen(bloom->vector));
      perror ("malloc");
      return -1;
    }
  else
    memset (bloom->vector, 0,
	    sizeof (char) * ((long long) (bloom->stat.elements / 8) + 1));
  /* generate a collection of random integers, to use later
   * when salting our keys before hashing them */
  //sketchy_randoms(&bloom->random_nums,hashes);
  //bloom->vector = "11111111";
  //printf("vector size-> %d\n",sizeof(bloom->vector));
  //memset(bloom->vector,0,sizeof(bloom->vector));
  return 0;
}

void
bloom_destroy (bloom * bloom)
{

  //free(bloom->random_nums.num);
  memset (bloom->vector, 0,
	  sizeof (char) * ((long long) (bloom->stat.elements / 8) + 1));
  free (bloom->vector);

}

int
bloom_check (bloom * bloom, char *str)
{
//printf("In bloom_check\n");
  return bloom_test (bloom, str, RO);
}

int
bloom_add (bloom * bloom, char *str)
{
  int ret;
  //printf("key--> %s\n",str);
  ret = bloom_test (bloom, str, SET);
  if (ret == 0)
    {
      bloom->inserts++;
    }
  return ret;
}

int
bloom_test (bloom * bloom, char *str, int mode)
{
  int i, hit;
  BIGNUM ret;
  //printf("In test\n");
  /* as many times as our ideal hash count dictates, salt our key
   * and hash it into the bit vector */
  hit = 1;
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

BIGNUM
bloom_hash (bloom * bloom, char *str, int i, int length)
{
  BIGNUM ret = 0;
  BIGNUM hash;

  ret = (BIGNUM) hash5 (str, seed[i], length) % (BIGNUM) bloom->stat.elements;

  return ret;
}

int
set (char *big, BIGNUM index)
{
  deref dr;

  finder (index, &dr);
  big[dr.index] += dr.spot;

  return 0;
}

int
test (char *big, BIGNUM index)
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

int
finder (BIGNUM index, deref * dr)
{

  dr->index = (BIGNUM) (index / 8);
  dr->spot = pow (2, (index % 8));
  //dr->index = (BIGNUM) (index >> 3);
  //dr->spot = pow (2, (index % 8));
  //dr->spot = 0x80;
  //dr->spot = dr->spot >> (index & 0x07);
  //dr->spot = pow(2,(index & 0x07));
  return 0;
}



int
sketchy_randoms (randoms * rands, int cnt)
{
  int i;
  srand (CONS);
  if ((rands->num = (int *) malloc (sizeof (int) * (cnt + 1))) == NULL)
    {
      perror ("malloc");
      errno = ERR_MALLOC;
      return -1;
    }
  for (i = 0; i < cnt; i++)
    {
      rands->num[i] =
	1 + (int) ((float) SALTRAND * rand () / (RAND_MAX + 1.0));
    }
  rands->cnt = cnt;
  return 0;
}

BIGNUM
report_capacity (bloom * bloom)
{
  return bloom->stat.capacity;
}

char *
to_bitstr (bloom * bm)
{
  int i, j;
  char *ptr;
  BIGNUM steps = (bm->stat.elements / 8) + 1;
  char *ret = (char *) malloc (sizeof (char *) * (bm->stat.elements + 2));

  strcpy (ret, "");
  for (i = 0, ptr = bm->vector; i < steps; i++, ptr++)
    {
      for (j = 0; j < 8; j++)
	{
	  if ((*ptr >> j) == 1)
	    {
	      strcat (ret, "1");
	    }
	  else
	    {
	      strcat (ret, "0");
	    }
	}
    }
  return ret;
}

int
save_bloom (char *filename, bloom * bl, char *prefix, char *target)
{
  char *position1,*position2;
  char *bloom_file = (char *) malloc (300 * sizeof (char));

  memset (bloom_file, 0, 300);

#ifdef DEBUG
  printf("target->%s\n",target);
#endif

  position1 = strrchr (filename, '/');
  position2 = strrchr (target+2, '/');

  if (prefix)
    strcat (bloom_file, prefix);
  else
    {
      if (position2)
      strncat (bloom_file,target+2,position2+1-(target+2));   
    
      if (position1)
	strncat (bloom_file, position1 + 1,
		 strrchr (filename, '.') - position1);
      else
	strncat (bloom_file, filename,
		 strrchr (filename, '.') - filename + 1);
      strcat (bloom_file, "bloom");
    }
  printf ("bloom name->%s\n", bloom_file);

  int fd, fq;

#ifdef __APPLE__
  fd = open(bloom_file, O_RDWR|O_CREAT, 0644);
#else // assume linux
  fd = open(bloom_file, O_RDWR|O_CREAT|O_LARGEFILE, 0644);
#endif
  if (fd < 0)
    {
      perror (bloom_file);
      return -1;
    }

  BIGNUM total_size =
    sizeof (bloom) + sizeof (char) * ((long long) (bl->stat.elements / 8) +
				      1) +
    sizeof (int) * (bl->stat.ideal_hashes + 1);

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

  fq = write (fd, bl, sizeof (bloom));

  total_size = (long long) (bl->stat.elements / 8) + 1;

  BIGNUM off = 0;
  while (total_size > TWOG)
    {
      fq = write (fd, bl->vector + off, sizeof (char) * TWOG);
      total_size -= TWOG;
      off += TWOG;
    }
  fq = write (fd, bl->vector + off, sizeof (char) * total_size);
  close (fd);

  memset (bl->vector, 0,
	  sizeof (char) * ((long long) (bl->stat.elements / 8) + 1));

  printf ("big file process OK\n");

  return 1;

}

int
load_bloom (char *filename, bloom * bl)
{
  int fd = 0;

  int x;

  printf ("bloom name->%s\n", filename);

#ifdef __APPLE__
  fd = open(filename, O_RDONLY, 0644); 
#else
  fd = open64(filename, O_RDONLY, 0644); 
#endif
  x = read (fd, bl, sizeof (bloom));

  bl->vector =
    (char *) malloc (sizeof (char) *
		     ((long long) (bl->stat.elements / 8) + 1));

  BIGNUM off = 0, total_size = ((long long) (bl->stat.elements / 8) + 1);

  while (total_size > TWOG)
    {
      x = read (fd, bl->vector + off, sizeof (char) * TWOG);

      total_size -= TWOG;

      off += TWOG;
    }

  x = read (fd, bl->vector + off, sizeof (char) * total_size);
  close (fd);

  printf ("successful bloom read...\n");

  close (fd);
}

void
write_result (char *filename, char *detail)
{
  int fd, x;

  fd = open (filename, O_CREAT | O_RDWR, S_IRWXU);

  x = write (fd, detail, strlen (detail));

  close (fd);
}

void
rev_trans (char *s)
{

  int i;
  int j;

  for (i = 0, j = strlen (s) - 1; i < j; ++i, --j)
    {
      char c = s[i];
      s[i] = s[j];
      s[j] = c;
    }

  i = 0;

  while (i < strlen (s))
    {
      switch (s[i])
	{
	case 'A':
	  s[0] = 'T';
	  break;
	case 'C':
	  s[0] = 'G';
	  break;
	case 'G':
	  s[0] = 'C';
	  break;
	case 'T':
	  s[0] = 'A';
	  break;
	case 'a':
	  s[0] = 't';
	  break;
	case 'c':
	  s[0] = 'g';
	  break;
	case 'g':
	  s[0] = 'c';
	  break;
	case 't':
	  s[0] = 'a';
	  break;
	}
      s++;
    }

}

/*
void instruction(){
printf("pod\n");



printf("SYNOPSIS\n");

printf("\n");

printf("ARGUMENTS\n");

printf("./script <>");

printf("The word size for K-mers used by the filter.\n");

printf("<bloomfilterlistfile>\n");

printf("File containing a list of Bloom filters built using I<bloombuild>\n");

printf("<queryfile>\n");

printf("Fasta file with reads to be filtered.\n");

printf("<outprefix>\n");

printf("Prefix for output filenames.\n");

printf("\n");

printf("DESCRIPTION\n");

printf("facs interrogates queries against filters and classifies queries onto genomes. The algorithm loops trough all queries for one filter at a time.\n");

printf("Results are written to three files\n");

exit(1);

}
*/
void
help ()
{
  printf
    ("##########################################################################\n");
  printf ("#NAME\n");
  printf ("#facs - Filter reads of DNA\n");
  printf ("# Contamination Removal tool\n");
  printf ("#\n");
  printf
    ("# extended from FACS (Fast and Accurate Classification of Sequences)\n");
  printf ("#\n");
  printf ("# facs.scilifelab.se\n");
  printf ("#\n");
  printf ("# Author Enze Liu\n");
  printf ("#\n");
  printf ("# enze.liu@scilifelab.se\n");
  printf
    ("##########################################################################\n");
  printf ("\n");
  printf ("DESCRIPTION\n");
  printf ("\n");
  printf
    ("    Bloom filters are a lightweight duplicate detection algorithm proposed\n");
  printf ("    by Burton Bloom\n");
  printf
    ("    (http://portal.acm.org/citation.cfm?id=362692&dl=ACM&coll=portal), with\n");
  printf
    ("    applications in stream data processing, among others. Bloom filters are\n");
  printf
    ("    a very cool thing. Where occasional false positives are acceptable,\n");
  printf
    ("    bloom filters give us the ability to detect duplicates in a fast and\n");
  printf ("    resource-friendly manner.\n");
  printf ("\n");
  printf
    ("    The allocation of memory for the bit vector is handled in the c layer,\n");
  printf
    ("    but perl's oo capability handles the garbage collection. when a\n");
  printf
    ("    Bloom::Faster object goes out of scope, the vector pointed to by the c\n");
  printf
    ("    structure will be free()d. to manually do this, the DESTROY builtin\n");
  printf ("    method can be called.\n");
  printf ("\n");
  printf
    ("    A bloom filter perl module is currently avaible on CPAN, but it is\n");
  printf
    ("    profoundly slow and cannot handle large vectors. This alternative uses a\n");
  printf
    ("    more efficient c library which can handle arbitrarily large vectors (up\n");
  printf ("    to the maximum size of a long long datatype (at least\n");
  printf ("    9223372036854775807, on supported systems ).\n");
  printf ("\n");
  printf
    ("    FACS is a novel algorithm, using Bloom filter to accurately and rapidly\n");
  printf
    ("    align sequences to a reference sequence. FACS was first optimized and\n");
  printf
    ("    validated using a synthetic metagenome dataset. An experimental metagenome\n");
  printf
    ("    dataset was then used to show that FACS is at least three times faster and\n");
  printf
    ("    more accurate than BLAT and SSAHA2 in classifying sequences when using\n");
  printf ("    references larger than 50Mbp.\n");
  printf ("\n");
}

void
build_help ()
{
  printf ("USAGE\n");
  printf
    ("##########################################################################\n");
  printf ("#   For Bloom build:\n");
  printf ("#   ./simple_check -m 1 -k 21 -e 0.005 -l file_list -p test/\n");
  printf ("#\n");
  printf ("#   Arguments:\n");
  printf ("#   -m mode (default 1)\n");
  printf ("#   1--> build one filter for each file in the list\n");
  printf ("#   2--> build one filter for every file in the list\n");
  printf ("#   -k k_mer size (default 1)\n");
  printf ("#   -e error rate (default 0.0005)\n");
  printf ("#   -l list containing all references files. One file per line\n");
  printf ("#   -p prefix (default same path as the script)\n");
  printf
    ("##########################################################################\n");
  exit (1);
}

void
check_help ()
{
  printf ("USAGE\n");
  printf
    ("##########################################################################\n");
  printf ("#   For Bloom check:\n");
  printf
    ("#   ./simple_check -m 1 -t 0.8 -s 1 -q test.fna -l bloom_list -p test/\n");
  printf ("#\n");
  printf ("#   Arguments:\n");
  printf ("#   -m mode (default 1)\n");
  printf ("#   1--> Bloom check\n");
  printf ("#   2--> Currently under evaluation\n");
  printf ("#   -t tolerant rate (default 0.8)\n");
  printf
    ("#   tolerant_rate is a threshold number between 0 and 1. Any read with a hit score more than tolerant_rate with be considered \n");
  printf ("#   as captured.\n");
  printf ("#   -s sampling rate (default 1)\n");
  printf
    ("#   sample_rate is a percentage. e.g. 0.2 means you only want to exam 20% reads out of total.\n");
  printf ("#   -q query file.\n");
  printf ("#   -l list containing all references files. One file per line\n");
  printf ("#   -p prefix (default same path as the script)\n");
//printf("#\n");   
//printf("#   *'1' is mode 1. For instance, when you use a ecoli filter and want to capture every contaminated read that caused by\n"); 
//printf("#   ecoli in the 'human.fna' query file, use mode 1. Mode 2 is currently under evaluation\n");
//printf("# \n");
//printf("#   *When you would like to remove every possible contamination in human.fna, you can use a human bloom filter and mode 2 \n");
//printf("#   to do that. It will capture everything that doesn't belong to human.\n");
//printf("#\n");
  printf
    ("##########################################################################\n");
  exit (1);
}

void
remove_help ()
{
  printf ("USAGE\n");
  printf
    ("##########################################################################\n");
  printf ("#   For Bloom remove:\n");
  printf
    ("#   ./simple_remove -m 1 -t 0.8 -q test.fna -l bloom_list -p test/\n");
  printf ("#\n");
  printf ("#   Arguments:\n");
  printf ("#   -m mode (default 1)\n");
  printf ("#   1--> Bloom remove\n");
  printf ("#   2--> Currently under evaluation\n");
  printf ("#   -t tolerant rate (default 0.8)\n");
  printf
    ("#   tolerant_rate is a threshold number between 0 and 1. Any read with a hit score more than tolerant_rate with be considered \n");
  printf ("#   as captured.\n");
  printf ("#   -q query file.\n");
  printf ("#   -l list containing all references files. One file per line\n");
  printf ("#   -p prefix (default same path as the script)\n");
  printf
    ("##########################################################################\n");
  exit (1);
}
