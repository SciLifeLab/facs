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


static char *
get_file_path (const char *path, const char *file_name)
{
  char *file_path = (char *) malloc (strlen (path) + 1 + strlen (file_name) + 1);
  strcpy (file_path, path);
  if (file_path[strlen (file_path) - 1] != '/')
    strcat (file_path, "/");
  strcat (file_path, file_name);
  return file_path;
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

  char *mimi = NULL;
  char *pos = NULL;
  int number = 0;

  if (list_user)
    {
      list_user = mmaping (list_user);
      while (list_user != '\0')
	{

	  mimi = (char *) malloc (200 * sizeof (char) + 1);
	  memset (mimi, 0, 200);

	  F_set *fset = NEW (F_set);
	  fset->filename = NULL;
	  fset->number = 0;

	  if ((pos = strchr (list_user, '\n')) != NULL)
	    strncat (mimi, list_user, pos - list_user);
	  else
	    break;

	  fset->filename = mimi;
	  int num = number;
	  fset->number = num;
	  fset->next = head->next;
	  head->next = fset;
	  head = head->next;
	  list_user = pos + 1;
	  number++;
	}
    }

  else if (is_file (file_user))
    {

      F_set *fset = NEW (F_set);
      fset->filename = file_user;
      fset->next = head->next;
      head->next = fset;
      head = head->next;
      head->next = NULL;
    }

  else if (is_dir (file_user))
    {
      dir = opendir (file_user);
      if (dir == NULL)
	{
	  perror ("Empty dir\n");
	  exit (-1);
	}
      while ((dir_info = readdir (dir)) != NULL)
	{
	  if (is_special_dir (dir_info->d_name))
	    continue;

	  if (!strstr (dir_info->d_name, ".bloom"))
	    continue;

	  F_set *fset = NEW (F_set);
	  fset->filename = get_file_path (file_user, dir_info->d_name);
	  fset->number = &number;
	  printf ("fset->%d\n", fset->number);
	  fset->next = head->next;
	  head->next = fset;
	  head = head->next;
	  number++;
	}
    }
  else
    {
      perror (file_user);
      exit (-1);
    }
  //free(file_path);
  //free(mimi);
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
