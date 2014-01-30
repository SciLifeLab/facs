#ifndef _BIGQUERY
#define _BIGQUERY

//#include "check.h"
#include "bloom.h"
#include <zlib.h>
extern BIGCAST CHUNKer(gzFile zip, BIGCAST offset, int chunk, char *data, char type);
extern BIGCAST CHUNKgz(gzFile zip, BIGCAST offset, int chunk, char *position, char *extra, char type);
extern int bq_main(int argc, char **argv);
extern char *query(char *query, char *bloom_filter, double tole_rate, double sampling_rate, char *list, char *target_path, char* report_fmt, char mode);
extern char *report(F_set * File_head, char* query, char* fmt, char* prefix, char* start_timestamp, double prob, int threads);
extern char *statistic_save (char *filename, char *prefix);
extern char *remove_data(char *data, char *start, char *end);
extern char *re_clean();
extern char *re_contam();
char *get_abs_path(char *filename);
extern char *move_start_point(char *filename);
extern char *strrstr(char *s, char *str);
extern void reset_string();
extern void init_string(int chunk);
extern void clean_list(Queue *head, Queue *tail);
extern void read_process(bloom *bl, Queue *info, Queue *tail, F_set *File_head, float sampling_rate, float tole_rate, char mode, char fmt_type);
#endif
