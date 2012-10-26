#include "bloom.h"
extern void evaluate (char *detail, char *filename, F_set *File_head);
extern void statistic_save (char *detail, char *filename, char* prefix);
extern void fasta_process (bloom * bl, Queue * info, Queue * tail, F_set *File_head, float sampling_rate, float tole_rate);
extern void fastq_process (bloom * bl, Queue * info, Queue * tail, F_set *File_head, float sampling_rate, float tole_rate);
extern int check_main(char *source, char *ref, float tole_rate, float sampling_rate, char *list, char *prefix, int help);
//extern int check(char *query, char *reference, char l, char *target_path, float sampling_rate, float tole_rate);
//extern int checky (char *query, char *reference, char l, char *target_path, float sampling_rate, float tole_rate);
