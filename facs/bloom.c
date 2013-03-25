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
#include "file_dir.h"

/*---------------------------*/
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>

int seed[20] =
  { 152501029, 152501717, 152503097, 152500171, 152500157, 152504837,
  10161313, 10371313, 10431313, 10501313, 10581313, 10611313, 10641313,
  10651313,
  10671313, 10731313, 10821313, 10881313, 10951313, 11001313
};

int
bloom_init (bloom * bloom, BIGNUM size, BIGNUM capacity, double error_rate,
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
  fprintf (stderr, "bloom_init(%lld,%d) => (%lld,%d) =>%f\n",
	   (BIGCAST) size, hashes, (BIGCAST) bloom->stat.elements,
	   bloom->stat.ideal_hashes, bloom->stat.e);
#endif

  if ((size > TOPLIMIT))
    {
      fprintf (stderr, "overflow2\n");
      return -2;
    }

  /* allocate our array of bytes.  where m is the size of our desired 
   * bit vector, we allocate m/8 + 1 bytes. */
  if ((bloom->vector = (char *) malloc (sizeof (char) *
					((long long) (bloom->stat.elements /
						      8) + 1))) == NULL)
    {
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

  memset (bloom->vector, 0,
	  sizeof (char) * ((long long) (bloom->stat.elements / 8) + 1));
  free (bloom->vector);
  bloom->vector = NULL;
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

  //dr->index = (BIGNUM) (index / 8);
  //dr->spot = 1<<(2, (index % 8));
  dr->index = (BIGNUM) (index >> 3);
  //dr->spot = pow (2, (index % 8));
  //dr->spot = 0x80;
  //dr->spot = dr->spot >> (index & 0x07);
  //dr->spot = pow(2,(index & 0x07));
  dr->spot = 1 << (index & 0x07);
  return 0;
}


BIGNUM
report_capacity (bloom * bloom)
{
  return bloom->stat.capacity;
}

char *
prefix_make (char *filename, char *prefix, char *target)
{
  char *position1 = strrchr (filename, '/');

  char *bloom_file = (char *) malloc (300 * sizeof (char));
  memset (bloom_file, 0, 300);
  if (is_dir (target))
    {
      strcat (bloom_file, target);
      if (position1 != NULL)
	strncat (bloom_file, position1, strrchr (position1, '.') - position1);
    }
  else if (target)
    {
      strcat (bloom_file, target);
    }
  else if (target != NULL && prefix != NULL)
    {
      if (position1 != NULL)
	strncat (bloom_file, position1, strrchr (position1, '.') - position1);
      else
	strncat (bloom_file, filename, strrchr (filename, '.') - filename);
      strcat (bloom_file, ".bloom");
      bloom_file++;
    }
  else
    {
      if (position1 != NULL)
	strncat (bloom_file, position1, strrchr (position1, '.') - position1);
      else
	strncat (bloom_file, filename, strrchr (filename, '.') - filename);
    }

  return bloom_file;
}

int
save_bloom (char *filename, bloom * bl, char *prefix, char *target)
{
  char *bloom_file = NULL;
  int fd;

  bloom_file = prefix_make (filename, prefix, target);

#ifdef DEBUG
  printf ("Bloom file to be written in: %s\n", bloom_file);
#endif

  if (prefix == NULL && target == NULL)
    {
      strcat (bloom_file, ".bloom");
      bloom_file++;
    }
  else if (is_dir (target))
    strcat (bloom_file, ".bloom");

#ifdef __APPLE__
  fd = open (bloom_file, O_RDWR | O_CREAT, PERMS);
#else // assume linux
  fd = open (bloom_file, O_RDWR | O_CREAT | O_LARGEFILE, PERMS);
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
  if (ftruncate (fd, total_size) < 0)
#else
  if (ftruncate64 (fd, total_size) < 0)
#endif
    {
      printf ("[%d]-ftruncate64 error: %s/n", errno, strerror (errno));
      close (fd);
      return 0;
    }

  if (write (fd, bl, sizeof (bloom)) < 0)
    {
      perror (" error writing bloom file ");
      exit (EXIT_FAILURE);
    }

  total_size = (long long) (bl->stat.elements / 8) + 1;

  BIGNUM off = 0;
  while (total_size > TWOG)
    {
      if (write (fd, bl->vector + off, sizeof (char) * TWOG) < 0)
	{
	  perror (" error writing bloom file ");
	  exit (EXIT_FAILURE);
	}
      total_size -= TWOG;
      off += TWOG;
    }
  if (write (fd, bl->vector + off, sizeof (char) * total_size) < 0)
    {
      perror (" error writing bloom file ");
      exit (EXIT_FAILURE);
    };
  close (fd);

  memset (bl->vector, 0,
	  sizeof (char) * ((long long) (bl->stat.elements / 8) + 1));

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
      exit (-1);
    }

  if (read (fd, bl, sizeof (bloom)) < 0)
    {
      perror ("Problem reading bloom filter");
    };

  bl->vector =
    (char *) malloc (sizeof (char) *
		     ((long long) (bl->stat.elements / 8) + 1));

  BIGNUM off = 0, total_size = ((long long) (bl->stat.elements / 8) + 1);

  while (total_size > TWOG)
    {
      ret = read (fd, bl->vector + off, sizeof (char) * TWOG);
      if (ret < 0)
	perror ("Problem reading bloom filter");
      total_size -= TWOG;
      off += TWOG;
    }

  ret = read (fd, bl->vector + off, sizeof (char) * total_size);

#ifdef DEBUG
  if (ret > 0)
    printf ("bloom filter read successfully\n");
  else
    ret = errno;
#endif

  close (fd);
  return ret;
}

void
write_result (char *filename, char *detail)
{
  int fd;

  fd = open (filename, O_CREAT | O_RDWR, S_IRWXU);
  if (write (fd, detail, strlen (detail)) < 0)
    {
      perror (" error writing result file ");
      exit (EXIT_FAILURE);
    }

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

char *
mmaping (char *source)
{
  struct stat statbuf;

  int fd = 0;
  char *sm = NULL;

  if ((fd = open (source, O_RDONLY | O_LARGEFILE)) < 0) {
      fprintf (stderr, "%s: %s\n", source, strerror (errno));
      exit (EXIT_FAILURE);
  }

  if (fstat (fd, &statbuf) < 0) {
      fprintf (stderr, "%s: %s\n", source, strerror (errno));
      exit (EXIT_FAILURE);
  } else if (statbuf.st_size == 0) {
      fprintf (stderr, "%s: %s\n", source, "File is empty");
      exit (-1);
  }

  sm = mmap (0, (BIGCAST) statbuf.st_size, PROT_READ,
	     MAP_SHARED | MAP_NORESERVE, fd, 0);

  if (MAP_FAILED == sm) {
      fprintf (stderr, "%s: %s\n", source, strerror (errno));
      exit (EXIT_FAILURE);
  }

  return sm;
}

char *
large_load (char *fifoname)
{
  int x = 0;
  char ch;
  FILE *fd;

  printf ("fifoname->%s\n", fifoname);
#ifdef __APPLE__
  fd = fopen (fifoname, "r");
#else
  fd = fopen64 (fifoname, "r");
#endif

  char *data = (char *) malloc ((TWOG / 2 + 1) * sizeof (char));

  data[TWOG / 2] = '\0';

  while ((ch = fgetc (fd)) != EOF)
    {
      data[x] = ch;
      x++;
    }

#ifdef DEBUG
  printf ("data length->%lld\n", (long long int) strlen (data));
#endif

  fclose (fd);
  return data;
}

BIGCAST
get_size (char *filename)
{
  struct stat buf;
  if (stat (filename, &buf) != -1)
    return buf.st_size;
  else
    return 0;
}
