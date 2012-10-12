#include <Python.h>
#include "bloom.h"
#include "file_dir.h"

static char module_docstring[] =
    "This module provides an interface for building and querying DRASS bloom filters";
static char bloom_docstring[] =
    "Builds a DRASS bloom filter and performs queries against it.";

/* Available functions */
static PyObject *drass_bloom_build(PyObject *self, PyObject *args, PyObject *argv);

static PyMethodDef module_methods[] = {
        {"build", drass_bloom_build, METH_VARARGS | METH_KEYWORDS, bloom_docstring},
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

   char *source, *bloom_filter;

   if (!PyArg_ParseTuple(args, "ss", &source, &bloom_filter))
       return NULL;
   
   build(source,21,5,1000, bloom_filter);

   return Py_BuildValue("");
}
