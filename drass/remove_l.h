#ifndef _REMOVE_CONTAM_L
#define _REMOVE_CONTAM_L

#include "bloom.h"
extern void fasta_process_ml (F_set * File_head, bloom * bl, Queue * info, Queue * tail, char *clean, char *contam, float tole_rate);
extern void fastq_process_ml (F_set * File_head, bloom * bl, Queue * info, Queue * tail, char *clean, char *contam, float tole_rate);
extern void save_result_ml (char *source, char *obj_file, char *data, char *data2, int flag, int type, char* prefix);
extern void all_save (F_set * File_head2, Queue * head2, Queue * tail, char *source, char *clean, char *clean2, char *contam, char *contam2, char *position, int type, char *prefix);
extern int count_read (char *dick, char *next, int type);
extern int remove_main_l(float tole_rate, char *source, char *ref, char *list, char *prefix, int help);

#endif
