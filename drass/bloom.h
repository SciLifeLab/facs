#ifndef _BLOOM
#define _BLOOM

#include <limits.h>

#define BIGNUM unsigned long long
#define BIGNUM_STR "unsigned long"
#define BIGCAST long long
#define TOPLIMIT LONG_MAX
#define PERMS 0644
#define HUN 1000
#define NEW(type) (type *) malloc(sizeof(type))

#define FILE_MODE (S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)
#define TWOG 2000000000
#define hashsize(n) ((BIGNUM)1<<(n))
#define hashmask(n) (hashsize(n) - 1)

#if !HAVE_SQRTL
#define sqrtl(val) ((long double)sqrt((double)val))
#endif

#define SALTRAND 99999999
#define CONS 1234567

#define HASH_FUNC hash5

typedef struct
{
	BIGNUM index;
        char spot;

} deref;

typedef BIGNUM (*hash_t)(char *str);

typedef struct 
{
	int *num;
	int cnt;
} randoms;

typedef struct bloomstat
{
	BIGNUM elements; /* size of array */
	int ideal_hashes; /* num hash functions */
	BIGNUM capacity; /* number of elements */
	double e; /* max error rate */
} bl_stat;

typedef struct
{
	char *vector;
	hash_t hash;
 	BIGNUM inserts;
  struct bloomstat stat;
  int k_mer;
  int dx;
  float mcf;
} bloom;

typedef struct info
{
     char *location;
     short *score;
     short *number;                 			
     struct info *next;          
} Queue;

typedef struct file_list
{
	char *filename;
        short number;
        BIGCAST reads_num;
        BIGCAST reads_contam;
	struct file_list *next;
} F_set;
/* these are modes to test_all() */
#define RO 0
#define SET 1
#define BVERBOSE 2

/* errors */

#define ERR_MALLOC 1
#define ERR_OVERFLOW 2
#define ERR_UNKNOWN 3

BIGNUM mkprime(BIGNUM startval);


extern int bloom_init(bloom *bloom, BIGNUM size, BIGNUM capacity,
                      double error_rate, int hashes, hash_t hash, int flags);

extern int bloom_check(bloom *bloom,char *str);
extern int bloom_add(bloom *bloom,char *str);
extern int bloom_test(bloom *bloom,char *str,int MODE);
extern void bloom_destroy(bloom *bloom);

extern int sketchy_randoms(randoms *rands,int cnt);
extern int finder (BIGNUM index,deref *dr);
extern int set(char *big,BIGNUM index);
extern int test (char *big, BIGNUM index);
extern BIGNUM bloom_hash(bloom *bloom,char *str, int i, int length);
extern int bloom_hash_old(bloom *bloom,char *str, int i);

extern BIGNUM find_close_prime(BIGNUM m);
extern int get_suggestion(struct bloomstat *stats, BIGNUM n,double e);
extern BIGCAST get_size (char *filename);
extern int kmer_suggestion (BIGCAST size);
extern float mco_suggestion (int k_mer);
extern int is_prime(BIGNUM m);
extern void get_rec (struct bloomstat *stat);
extern BIGNUM report_capacity(bloom *bloom);

extern void write_result (char *filename, char *detail);
extern void build_help(void);
extern void check_help(void);
extern void remove_help(void);
extern void remove_l_help(void);
extern int save_bloom (char *filename, bloom *bl, char *prefix, char *target);
extern int load_bloom (char *filename, bloom *bl);
extern void rev_trans (char *s);
//extern void cat_print(char *merge, char *remove);

//extern void kmer_suggestion (BIGCAST size);
//extern float mcf_suggestion (int k_mer);
extern char *large_load (char *fifoname);
extern char *mmaping (char *source);
extern char *prefix_make (char *filename, char *prefix, char *target);
#endif
