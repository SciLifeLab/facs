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
   const *flags;
   char *source, *position, *prefix;
   struct stat statbuf;
   bloom *bl_2;
   bl_2 = NEW (bloom);
   position = mmaping(source);

/* initializes argc & argv, needs ParseTuple for each parm? */
   init(argc, argv);
   if (!PyArg_ParseTuple(args, "s", &flags))
       return NULL;

   if (*position == '>')
       capacity = strlen (position);
   else
       capacity = strlen (position) / 2;

   init_bloom (bl_2);
   fasta_add(bl_2, position);

   save_bloom(source, bl_2, prefix, flags);

   munmap (position, statbuf.st_size);
}
