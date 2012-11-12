#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include "file_dir.h"
#include "bloom.h"


void
get_file_path (const char *path, const char *file_name, char *file_path)
{
  strcpy (file_path, path);
  if (file_path[strlen (path) - 1] != '/')
    strcat (file_path, "/");
  strcat (file_path, file_name);
}

int
is_dir (const char *path)
{
  struct stat statbuf;
  if (lstat (path, &statbuf) == 0)
    {
      return S_ISDIR (statbuf.st_mode) != 0;
    }
  return 0;
}

int
is_file (const char *path)
{
  struct stat statbuf;
  if (lstat (path, &statbuf) == 0)
    return S_ISREG (statbuf.st_mode) != 0;
  return 0;
}

int
is_special_dir (const char *path)
{
  return strcmp (path, ".") == 0 || strcmp (path, "..") == 0;
}

F_set *
make_list (char *file_user, char *list_user)
{
  struct DIR *dir;
  struct dirent *dir_info;

  F_set *head = NEW (F_set);
  F_set *head1 = head;
  int number = 0;

  if (list_user)
    {
  
      list_user = mmaping (list_user);
      char *pos;


      while (list_user != '\0')
	{

	  char *mimi = (char *) malloc (300 * sizeof (char) + 1);
	  memset (mimi, 0, 300);

	  F_set *crap = NEW (F_set);

	  if (pos = strchr (list_user, '\n'))
	    strncat (mimi, list_user, pos - list_user);
	  else
	    break;

	  crap->filename = mimi;
	  crap->number = &number;
	  crap->next = head->next;
	  head->next = crap;
	  head = head->next;
	  list_user = pos + 1;
	  number++;
	}
    }
 
  else if (is_file (file_user))
    {
   
      F_set *crap = NEW (F_set);
      crap->filename = file_user;
      crap->next = head->next;
      head->next = crap;
      head = head->next;
      head->next = NULL;
    }

  else if (is_dir (file_user))
    {


      if ((dir = opendir (file_user)) == NULL)
	{
	  perror ("Empty dir\n");
	  exit (-1);
	}
      while ((dir_info = readdir (dir)) != NULL)
	{
	  char *file_path = (char *) malloc (300 * sizeof (char));
	  memset (file_path, 0, 300);
	  get_file_path (file_user, dir_info->d_name, file_path);

	  if (is_special_dir (dir_info->d_name))
	    continue;

	  if (!strstr (dir_info->d_name, ".bloom"))
	    continue;
	  printf ("file_path->%s\n", file_path);
	  F_set *crap = NEW (F_set);
	  crap->filename = file_path;
	  crap->number = &number;
	  printf ("crap->%d\n", crap->number);
	  crap->next = head->next;
	  head->next = crap;
	  head = head->next;
	  number++;
	}
    }
  head1->next->reads_num = 0;
  head1->next->reads_contam = 0;
  return head1->next;

}

/*
void main()
{
//char *list = "test_l";
char *list;
char *single_file = "single_file";
char *dir = "/home/apple/xt/dir";

F_set *head;

head = make_list(single_file,list);

head = head->next;

while (head)
      {
      printf("position->%s\n",head->filename);
      //printf("next->%d\n",head->next);
      head = head->next;
      }

return;
}
*/
