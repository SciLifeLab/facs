#ifndef _BIGQUERY
#define _BIGQUERY

#include "check.h"
#include "bloom.h"
#include <zlib.h>
extern char *bac_2_n (char *filename);
extern char *strrstr (char *s, char *str);
//extern BIGCAST get_size (char *filename);
extern void clean_list (Queue * head, Queue * tail);
<<<<<<< HEAD
extern BIGCAST CHUNKer (gzFile zip, BIGCAST offset, int chunk, char *data, int type);
extern BIGCAST CHUNKgz (gzFile zip, BIGCAST offset, int chunk, char *position, char *extra, int type);
extern int bq_main (int argc, char **argv);
extern int query (char *query, char *bloom_filter, double tole_rate, double sampling_rate, char *list, char *target_path, char* report_fmt, char mode);
=======
extern BIGCAST CHUNKer (gzFile zip, BIGCAST offset, int chunk, char *data,
			int type);
extern BIGCAST CHUNKgz (gzFile zip, BIGCAST offset, int chunk, char *position,
			char *extra, int type);
extern char* bq_main (int argc, char **argv);
extern char* query (char *query, char *bloom_filter, double tole_rate,
                  double sampling_rate, char *list, char *target_path, char* report_fmt);
>>>>>>> 026539af0185236b7ad9f514bb5f9b0eaba67fd8
#endif
