#define NEW(type) (type *) malloc(sizeof(type))

typedef struct file_list
{
	char *filename;
	struct file_list *next;
    short number;
    short reads_num;
    short reads_contam;
} F_set;

extern void get_file_path(const char *path, const char *file_name,  char *file_path);
extern int is_dir(const char *path);
extern int is_file(const char *path);
extern int is_special_dir(const char *path);
extern F_set *make_list(char *file_user, char *list_user);
