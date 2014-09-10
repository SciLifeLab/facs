#include <zlib.h>
#include "bloom.h"
static int mpicheck_usage (void);
int gather (F_set *File_head, int total_proc, int proc_num);
BIGCAST struc_init (char *filename, int proc_num, int total_proc);
BIGCAST gz_mpi (gzFile zip, BIGCAST offset, BIGCAST left, char *data, char type);
