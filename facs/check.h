#ifndef _QUERY
#define _QUERY

#include "bloom.h"
extern void read_process(bloom * bl, Queue * info, Queue * tail, F_set * File_head, float sampling_rate, float tole_rate, char mode, char fmt_type);
extern char *report(F_set * File_head, char* query, char* fmt, char* prefix, int k_mer);
extern char *statistic_save (char *filename, char *prefix);
extern char *remove_data(char *data, char *start, char *end);
extern char *re_clean();
extern char *re_contam();
extern void reset_string();
//extern void init_string(int chunk);
#endif
