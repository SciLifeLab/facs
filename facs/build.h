#ifndef _BUILD
#define _BUILD

#include "bloom.h"
extern void ref_add (bloom * bl, char *position);
extern void fastq_add (bloom * bl, char *position);
extern void fasta_add (bloom * bl, char *position);
extern char *fasta_data (bloom * bl_2, char *data);
extern void init_bloom (bloom * bl, BIGNUM capacity, float error_rate,
			int k_mer, char *filename);
extern int build (char *ref_name, char *target_path, int k_mer,
		  double error_rate, char *prefix);
extern int build_main (int argc, char **argv);
#endif
