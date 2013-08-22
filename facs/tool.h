#ifndef _TOOL
#define _TOOL

#include "bloom.h"
int fq_read_length (char *data);
extern int dx_add (int k_mer);
extern int get_parainfo (char *full, Queue *head);
extern char *fastq_relocate (char *data, int offset, int length);
extern char *jump (char *target, char type, float sampling_rate);
extern char *get_right_sp (char *start_point ,char type);
extern char *check_fmt (Queue *info, Queue *tail, char *start_point, char type);
int fastq_full_check (bloom * bl, char *p, int length, float tole_rate, F_set *File_head);
int fasta_full_check (bloom * bl, char *begin, char *next, char model, float tole_rate, F_set *File_head);
extern int fastq_read_check (char *begin, int length, char mode, bloom * bl, float tole_rate, F_set *File_head);
extern int fasta_read_check (char *begin, char *next, char model, bloom * bl, float tole_rate, F_set *File_head);
void isodate(char* buf);

#endif
