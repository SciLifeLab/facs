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


double get_mu (BIGNUM num_hit, double prob)
{
       return ((double)num_hit)*prob;
}

double get_sigma (BIGNUM num_hit, double prob)
{
       return (double)num_hit*prob*(1-prob);
}

/*
float get_probability (BIGCAST hits, BIGCAST total, int k_mer)
{
double times = (double)total/(100*MB);
double prob = 0;
int rand_hit = 0;
//((k_mer/3)==0)&&(k_mer=k_mer,1)||(k_mer-(k_mer%3));
switch (k_mer)
{
  case 9:
  rand_hit = 10;
  break;
  case 12:
  rand_hit = 20;
  break;
  case 15:
  rand_hit = 40;
  break;
  case 18:
  rand_hit = 80;
  break;
  case 21:
  rand_hit = 100;
  break;
  default:
  printf ("cant handle this k_mer so far\n");
  exit(-1);
}
prob = (double)hits/times;
if (prob<rand_hit)
    return 0;
else
    return hits/total; 
}
*/


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
  if (size<10*MB)
     {
      //bl->k_mer = 15;
      //bl->mcf = 0.4;
      return 15;
     }
  else if (size<20*MB)
     {
      //bl->k_mer = 17;
      //bl->mcf = 0.4;
      return 16;
     }
  else if (size<50*MB)
     {
      //bl->k_mer = 17;
      //bl->mcf = 0.4;
      return 17;
     }
  else if (size<100*MB)
     {
      //bl->k_mer = 18;
      //bl->mcf = 0.3;
      return 18;
     }
  else if (size<500*MB)
     {
      return 19;
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
switch (k_mer)
  {
    case 15:
    	return 0.4;
    case 16:
        return 0.3;
    case 17:
        return 0.3;
    case 18:
        return 0.3; 
    case 19:
        return 0.4;
    case 20:
        return 0.3;
    default:
        return 0.4;
  }
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
