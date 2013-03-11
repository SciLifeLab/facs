#ifndef _FILE_DIR
#define _FILE_DIR

#include "bloom.h"
extern int is_dir (const char *path);
extern int is_file (const char *path);
extern int is_special_dir (const char *path);
extern F_set *make_list (char *file_user, char *list_user);
extern void get_file_path (const char *path, const char *file_name,
			    char *file_path);

#endif	/*  */
