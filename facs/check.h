#ifndef _QUERY
#define _QUERY

#include "bloom.h"
extern void report(char *detail, char *filename, F_set * File_head, char* query, char* fmt, char* prefix);
extern void fasta_process (bloom * bl, Queue * info, Queue * tail, F_set * File_head, float sampling_rate, float tole_rate);
extern void fastq_process (bloom * bl, Queue * info, Queue * tail, F_set * File_head, float sampling_rate, float tole_rate, char mode);
extern char *statistic_save (char *filename, char *prefix);
extern char *remove_data(char *data, char *start, char *end);
extern char *re_clean();
extern char *re_contam();
extern void init_string(int chunk);
#endif
