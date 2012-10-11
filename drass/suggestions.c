#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bloom.h"

/* dmitriy ryaboy */

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
