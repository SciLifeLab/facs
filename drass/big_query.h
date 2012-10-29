#include "bloom.h"
extern char *strrstr(char *s, char *str);
extern BIGCAST get_size (char *filename);
extern BIGCAST CHUNKer(gzFile zip,BIGCAST offset,int chunk,char *data,int type);
extern BIGCAST CHUNKgz(gzFile zip, BIGCAST offset,int chunk,char *position,char *extra,int type);
extern int bq_main(char *source, char *ref, float tole_rate, float sampling_rate, char *list, char *prefix, int help);
