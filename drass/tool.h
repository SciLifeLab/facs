extern int fastq_full_check (bloom * bl, char *p, int distance,  char *model, float tole_rate);
extern int fasta_full_check (bloom * bl, char *begin, char *next, char *model, float tole_rate);
extern int fastq_read_check (char *begin, int length, char *model, bloom * bl, float tole_rate);
extern int fasta_read_check (char *begin, char *next, char *model, bloom * bl, float tole_rate);
