#include <Python.h>
#include "bloom.h"
#define NEW(type) (type *) malloc(sizeof(type))

static char module_docstring[] =
    "This module provides an interface for building and querying DRASS bloom filters";
static char bloom_docstring[] =
    "Builds a DRASS bloom filter and performs queries against it.";

/* Available functions */
static PyObject *drass_bloom_build(PyObject *self, PyObject *args, PyObject *argv);

static PyMethodDef module_methods[] = {
        {"bloom_build", drass_bloom_build, METH_VARARGS | METH_KEYWORDS, bloom_docstring},
/*        {"bloom_query", drass_bloom_query, METH_VARARGS | METH_KEYWORDS, bloom_docstring}, */
        {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initdrass(void)
{
    PyObject *m = Py_InitModule3("drass", module_methods, module_docstring);
    if (m == NULL)
        return;
}

static PyObject *drass_bloom_build(PyObject *self, PyObject *args, PyObject *argv)
{
   BIGNUM capacity;
   char* bloom_filter;
   char* mmaping (char* bloom);
   char *source, *position, *prefix;
   struct stat statbuf;
/* int argc; */
   char **argv_drass;
   bloom *bl_2;
   bl_2 = NEW(bloom);


   if (!PyArg_ParseTuple(args, "ss", &source, &bloom_filter))
       return NULL;

   printf("From drass bloom: %s\n", bloom_filter);
   position = mmaping(source);
   printf("From drass: %s\n", source);
   printf("From drass after mmaping: %s\n", source);

/* initializes argc & argv, needs ParseTuple for each parm? */
/*   init(argc, argv_drass); */

   if (*position == '>')
       capacity = strlen (position);
   else
       capacity = strlen (position) / 2;

   printf("What da fuck!? Capacity is: %d\n", capacity);
   init_bloom(bl_2);
   fasta_add(bl_2, position);

   save_bloom(source, bl_2, prefix, &bloom_filter);

   munmap (position, statbuf.st_size);

   return NULL;
}
