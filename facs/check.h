#ifndef _QUERY
#define _QUERY

#include "bloom.h"
extern void evaluate (char *detail, char *filename, F_set *File_head);
extern void statistic_save (char *detail, char *filename, char* prefix);
void fasta_process (bloom * bl, Queue * info, Queue * tail, F_set *File_head, float sampling_rate, float tole_rate);

#endif
