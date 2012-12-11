#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bloom.h"
/*------------------------------*/
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>
#define MB 1048576


/* dmitriy ryaboy */

BIGCAST get_size (char *filename)
{
BIGCAST tim;
struct stat statbuf;
if ((tim=open(filename, O_RDONLY))<0)
    {
     printf("open file error...\n");
     exit(-1);
    }
fstat (tim, &statbuf);
return statbuf.st_size;
}

int kmer_suggestion (BIGCAST size)
{
  if (size<1*MB)
     {
      //bl->k_mer = 12;
      //bl->mcf = 0.3;
      return 12;
     }
  else if (size<20*MB)
     {
      //bl->k_mer = 15;
      //bl->mcf = 0.4;
      return 15;
     }
  else if (size<50*MB)
     {
      //bl->k_mer = 17;
      //bl->mcf = 0.4;
      return 17;
     }
  else if (size<200*MB)
     {
      //bl->k_mer = 18;
      //bl->mcf = 0.3;
      return 18;
     }
  else 
     {
      //bl->k_mer = 20;
      //bl->mcf = 0.3;
      return 20;
     }
}

float mco_suggestion (int k_mer)
{
  if (k_mer<15)
      return 0.3;
  if (k_mer<18)
      return 0.4;
  else
      return 0.3;
}

int
get_suggestion (struct bloomstat *stats, BIGNUM n, double e)
{
  stats->capacity = n;
  stats->e = e;
  get_rec (stats);

  return 0;
}

BIGNUM
find_close_prime (BIGNUM m)
{
  if ((m % 2) == 0)
    m += 1;

  while (!is_prime (m))
    {
      m += 2;
    }
  return m;
}

/*
Given the desired capacity and error rate, calculate the appropriate values
for number of hash functions and size of array
*/
void
get_rec (struct bloomstat *stat)
{
  /* assuming perfect number of cells, k directly depends on e */
  stat->ideal_hashes = (int) log (stat->e) / log (0.5);
  stat->elements =
    find_close_prime ((BIGNUM) 13 * stat->capacity *
		      (BIGNUM) stat->ideal_hashes / (BIGNUM) 9);
  /*
     recalculate k with the actual m, not the ideal 
     wouldn't need to if it wasn't prime, but that causes problems
     for hash algs
   */
  stat->ideal_hashes = 9 * stat->elements / (13 * stat->capacity);
}

int
is_prime (BIGNUM m)
{
  BIGNUM a = (BIGNUM) sqrtl ((long double) m);
  BIGNUM currval = 3;
  if (m % 2 == 0 && m != 2)
    return 0;
  while (m % currval != 0 && currval < a)
    {
      if (m % (currval + 2) == 0)
	return 0;
      if (m % (currval + 4) == 0)
	return 0;
      currval += 8;
    }
  return (int) m % currval;
}
