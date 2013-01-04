#include <Python.h>
#include "bloom.h"
#include "file_dir.h"
#include "tool.h"
#include "build.h"
#include "check.h"
#include "big_query.h"

static char module_docstring[] =
    "This module provides an interface for building and querying DRASS bloom filters";
static char bloom_docstring[] =
    "Builds a DRASS bloom filter and performs queries against it.";

/* Available functions */
static PyObject *drass_bloom_build(PyObject *self, PyObject *args);
static PyObject *drass_bloom_query(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
        {"build", drass_bloom_build, METH_VARARGS | METH_KEYWORDS, bloom_docstring},
        {"query", drass_bloom_query, METH_VARARGS | METH_KEYWORDS, bloom_docstring},
        {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initdrass(void)
{
    PyObject *m = Py_InitModule3("drass", module_methods, module_docstring);
    if (m == NULL)
        return;
}


static PyObject *drass_bloom_query(PyObject *self, PyObject *args)
{
   double sampling_rate=1;
   double tole_rate=0;
   char *qry, *bloom;
   int ret;

   if (!PyArg_ParseTuple(args, "ss|dd", &qry, &bloom, &tole_rate, &sampling_rate))
       return NULL;
   ret = query(qry, bloom, tole_rate, sampling_rate, NULL, NULL);

   return Py_BuildValue("i", ret);
}

static PyObject *drass_bloom_build(PyObject *self, PyObject *args)
{
   char *source, *bloom_filter, *prefix;
   int ret;

   //DRASS operational defaults
   int k_mer=0;
   double error_rate=0.0005;

   if (!PyArg_ParseTuple(args, "ss|ids", &source, &bloom_filter, &k_mer, &error_rate, &prefix))
       return NULL;
   
   ret = build(source, bloom_filter, k_mer, error_rate, prefix);

   return Py_BuildValue("i", ret);
}
