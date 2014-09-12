#include <Python.h>
#include "bloom.h"
#include "file_dir.h"
#include "tool.h"
#include "build.h"
#include "query.h"
#include "remove.h"

static char module_docstring[] =
"This module provides an interface for building and querying FACS bloom \
filters";

static char build_docstring[] = "Builds a FACS bloom filter and performs queries \
against it.\nUsage: facs.build(\"ecoli.fasta\", \"ecoli.bloom\")";
static char query_docstring[] = "Performs a query against a FACS bloom filter.\
\n Usage: facs.query(\"contaminated_sample.fastq.gz\", \"ecoli.bloom\")";
static char remove_docstring[] = "Removes contaminated sequences from a \
FASTA/FASTQ file using a FACS bloom filter. \
\n Usage: facs.remove(\"contaminated_sample.fastq.gz\", \"ecoli.bloom\")";

/* Available functions */
static PyObject *facs_bloom_build (PyObject * self, PyObject * args);
static PyObject *facs_bloom_query (PyObject * self, PyObject * args);
static PyObject *facs_bloom_remove (PyObject * self, PyObject * args);

static PyMethodDef module_methods[] = {
  {"build", facs_bloom_build, METH_VARARGS | METH_KEYWORDS, build_docstring},
  {"query", facs_bloom_query, METH_VARARGS | METH_KEYWORDS, query_docstring},
  {"remove", facs_bloom_remove, METH_VARARGS | METH_KEYWORDS, remove_docstring},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_facs (void)
{
  PyObject *m = Py_InitModule3 ("_facs", module_methods, module_docstring);
  if (m == NULL)
    return;
}


static PyObject *
facs_bloom_query(PyObject * self, PyObject * args)
{
  double sampling_rate = 1;
  double tole_rate = 0;
  char *qry, *bloom;
  char* report_fmt = "json";
  char* ret;
  if (!PyArg_ParseTuple(args, "ss|dds", &qry, &bloom, &tole_rate, &sampling_rate, report_fmt))
    return NULL;
  ret = query(qry, bloom, tole_rate, sampling_rate, NULL, NULL, report_fmt,'c');

  printf("%s\n", ret);

  return Py_BuildValue("s", ret);
}

static PyObject *
facs_bloom_build(PyObject * self, PyObject * args)
{
  char *source, *bloom_filter, *prefix;
  int ret;

  //FACS operational defaults
  int k_mer = 0;
  double error_rate = 0.005;

  if (!PyArg_ParseTuple(args, "ss|ids", &source, &bloom_filter,
                        &k_mer, &error_rate, &prefix))
    return NULL;

  ret = build(source, bloom_filter, k_mer, error_rate, prefix);

  return Py_BuildValue ("i", ret);
}


static PyObject *
facs_bloom_remove(PyObject * self, PyObject * args)
{
  double tole_rate = 0;
  char *src, *ref, *list, *prefix;
  char *report_fmt = "json";
  char *ret;

  if (!PyArg_ParseTuple
      (args, "ss|ssd", &src, &ref, &list, &prefix, &tole_rate))
    return NULL;

  //ret = remove_reads(src, ref, NULL, NULL, tole_rate);
  ret = query(src, ref, tole_rate, 1.000,  NULL, NULL, report_fmt , 'r');
  printf("%s\n",ret);
  return Py_BuildValue ("i", ret);
}
