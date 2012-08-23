#ifndef _BLOOM
#define _BLOOM


#include <limits.h>

#define BIGNUM unsigned long long
#define BIGNUM_STR "unsigned long"
#define BIGCAST long long
#define TOPLIMIT LONG_MAX

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

struct bloomstat
{
	BIGNUM elements; /* size of array */
	int ideal_hashes; /* num hash functions */
	BIGNUM capacity; /* number of elements */
	float e; /* max error rate */
} ;

typedef struct
{
	char *vector;
	hash_t hash;
 	BIGNUM inserts;
  //	int ideal_hashes;
        struct bloomstat stat;
	//randoms random_nums;
        int k_mer;
} bloom;

//int verbose;

typedef struct info
{
     int number;
     char *location;                     			
     struct info *next;          
} Queue;

/* these are modes to test_all() */
#define RO 0
#define SET 1
#define BVERBOSE 2

/* errors */

#define ERR_MALLOC 1
#define ERR_OVERFLOW 2
#define ERR_UNKNOWN 3

BIGNUM mkprime(BIGNUM startval);


extern int bloom_init(bloom *bloom,BIGNUM size,BIGNUM capacity, float error_rate,int hashes,hash_t hash,int flags);
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
extern int is_prime(BIGNUM m);
extern void get_rec (struct bloomstat *stat);
extern BIGNUM report_capacity(bloom *bloom);

extern void write_result (char *filename, char *detail);
extern void help();
extern void build_help();
extern void check_help();
extern void remove_help();
extern int save_bloom (char *filename, bloom *bl, char *prefix, char *target);
extern int load_bloom (char *filename, bloom *bl);
extern void rev_trans (char *s);
//extern void cat_print(char *merge, char *remove);
#endif
