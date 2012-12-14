#ifndef _HASHES
#define _HASHES

extern BIGNUM hash1 (char *str, int seed);
extern BIGNUM hash2 (char *str, int seed);
extern BIGNUM hash3 (char *str, int seed);
extern BIGNUM hash4 (char *str, int seed);
extern BIGNUM hash5 (char *str, int seed, int length);

#endif
