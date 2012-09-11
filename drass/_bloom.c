#include <Python.h>

static char module_docstring[] =
    "This module provides an interface for building and querying DRASS bloom filters";
static char bloom_docstring[] =
    "Builds a DRASS bloom filter and performs queries against it.";

static PyMethodDef module_methods[] = {
        {"bloom_build", drass_bloom_build, METH_VARARGS | METH_KEYWORDS, bloom_docstring},
        {"bloom_query", drass_bloom_query, METH_VARARGS | METH_KEYWORDS, bloom_docstring},
        {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_bloom(void)
{
    PyObject *m = Py_InitModule3("_bloom", module_methods, module_docstring);
    if (m == NULL)
        return;
}

static PyObject *drass_bloom_build(PyObject *self, PyObject *args);
static PyObject *drass_bloom_query(PyObject *self, PyObject *args);
