#ifndef _REMOVE_CONTAM
#define _REMOVE_CONTAM

#include "bloom.h"
extern void fasta_process_m (bloom * bl, Queue * info, Queue * tail,
			     float tole_rate, F_set * File_head, int type);
extern void fastq_process_m (bloom * bl, Queue * info, Queue * tail,
			     float tole_rate, F_set * File_head, int type);
extern int remove_main (int argc, char **argv);
extern int remove_reads (char *source, char *ref, char *list, char *prefix,
			 float tole_rate);
extern void save_result (char *source, char *obj_file, char type, char *prefix,
			 char *clean2, char *contam2);

#endif
