#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bloom.h"
#include "prob.h"
/*------------------------------*/
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>
#include "prob.h"
#define MB 1048576
/*------------------------------*/
/*
get_suggestion, find_close_prime, get_rec and is_prime functions are invoked from bloom::faster
*/
/*------------------------------*/
/*suggestion for probability of random hits based on different k-mer in a empiric way*/
double prob_suggestion (int k_mer)
{
double prob = 0;
if (k_mer<=12)
	prob = 0.51038;
else if (k_mer<=15)
	prob = 0.05569;
else if (k_mer<=18)
	prob = 0.00636;
else 
	prob = 0.001057;
return prob;
}
/*suggestion for kmer,suggest different k_mer number based on the size of the sample*/
int kmer_suggestion (BIGCAST size)
{
  if (size < 10 * MB)
    {
      return 15;
    }
  else if (size < 20 * MB)
    {
      return 16;
    }
  else if (size < 50 * MB)
    {
      return 17;
    }
  else if (size < 100 * MB)
    {
      return 18;
    }
  else if (size < 500 * MB)
    {
      return 19;
    }
  else
    {
      return 20;
    }
}
/*suggestion for match cutoff (threshold) based on the k-mer number*/
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
/*suggestion for size of bloom filter, estimated by the number of elements in the query*/
int get_suggestion (struct bloomstat *stats, BIGNUM n, double e)
{
  stats->capacity = n;
  stats->e = e;
  get_rec (stats);

  return 0;
}
/*find the closest prme number for number of elements*/
BIGNUM find_close_prime (BIGNUM m)
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
void get_rec (struct bloomstat *stat)
{
  /* assuming perfect number of cells, k directly depends on e */
  stat->ideal_hashes = (int) log (stat->e) / log (0.5);
  stat->elements = find_close_prime((BIGNUM) 13 * stat->capacity *(BIGNUM) stat->ideal_hashes / (BIGNUM) 9);
  /*
     recalculate k with the actual m, not the ideal 
     wouldn't need to if it wasn't prime, but that causes problems
     for hash algs
   */
  stat->ideal_hashes = 9 * stat->elements / (13 * stat->capacity);
}

int is_prime (BIGNUM m)
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
