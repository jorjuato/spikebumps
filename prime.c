/* Generated by Pyrex 0.9.4.1 on Wed May 16 04:07:42 2007 */

#include "Python.h"
#include "structmember.h"
#ifndef PY_LONG_LONG
  #define PY_LONG_LONG LONG_LONG
#endif
#ifdef __cplusplus
#define __PYX_EXTERN_C extern "C"
#else
#define __PYX_EXTERN_C extern
#endif
__PYX_EXTERN_C double pow(double, double);


typedef struct {PyObject **p; char *s;} __Pyx_InternTabEntry; /*proto*/
typedef struct {PyObject **p; char *s; long n;} __Pyx_StringTabEntry; /*proto*/
static PyObject *__Pyx_UnpackItem(PyObject *, int); /*proto*/
static int __Pyx_EndUnpack(PyObject *, int); /*proto*/
static int __Pyx_PrintItem(PyObject *); /*proto*/
static int __Pyx_PrintNewline(void); /*proto*/
static void __Pyx_Raise(PyObject *type, PyObject *value, PyObject *tb); /*proto*/
static void __Pyx_ReRaise(void); /*proto*/
static PyObject *__Pyx_Import(PyObject *name, PyObject *from_list); /*proto*/
static PyObject *__Pyx_GetExcValue(void); /*proto*/
static int __Pyx_ArgTypeTest(PyObject *obj, PyTypeObject *type, int none_allowed, char *name); /*proto*/
static int __Pyx_TypeTest(PyObject *obj, PyTypeObject *type); /*proto*/
static int __Pyx_GetStarArgs(PyObject **args, PyObject **kwds, char *kwd_list[], int nargs, PyObject **args2, PyObject **kwds2); /*proto*/
static void __Pyx_WriteUnraisable(char *name); /*proto*/
static void __Pyx_AddTraceback(char *funcname); /*proto*/
static PyTypeObject *__Pyx_ImportType(char *module_name, char *class_name, long size);  /*proto*/
static int __Pyx_SetVtable(PyObject *dict, void *vtable); /*proto*/
static int __Pyx_GetVtable(PyObject *dict, void *vtabptr); /*proto*/
static PyObject *__Pyx_CreateClass(PyObject *bases, PyObject *dict, PyObject *name, char *modname); /*proto*/
static int __Pyx_InternStrings(__Pyx_InternTabEntry *t); /*proto*/
static int __Pyx_InitStrings(__Pyx_StringTabEntry *t); /*proto*/
static PyObject *__Pyx_GetName(PyObject *dict, PyObject *name); /*proto*/

static PyObject *__pyx_m;
static PyObject *__pyx_b;
static int __pyx_lineno;
static char *__pyx_filename;
static char **__pyx_f;

/* Declarations from prime */


/* Implementation of prime */

static PyObject *__pyx_n_primes;

static PyObject *__pyx_n_append;

static PyObject *__pyx_f_5prime_primes(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); /*proto*/
static PyObject *__pyx_f_5prime_primes(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  int __pyx_v_kmax;
  int __pyx_v_n;
  int __pyx_v_k;
  int __pyx_v_i;
  int (__pyx_v_p[1000]);
  PyObject *__pyx_v_result;
  PyObject *__pyx_r;
  PyObject *__pyx_1 = 0;
  int __pyx_2;
  PyObject *__pyx_3 = 0;
  PyObject *__pyx_4 = 0;
  static char *__pyx_argnames[] = {"kmax",0};
  if (!PyArg_ParseTupleAndKeywords(__pyx_args, __pyx_kwds, "i", __pyx_argnames, &__pyx_v_kmax)) return 0;
  __pyx_v_result = Py_None; Py_INCREF(Py_None);

  /* "/root/prime.pyx":8 */
  __pyx_1 = PyList_New(0); if (!__pyx_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 8; goto __pyx_L1;}
  Py_DECREF(__pyx_v_result);
  __pyx_v_result = __pyx_1;
  __pyx_1 = 0;

  /* "/root/prime.pyx":9 */
  __pyx_2 = (__pyx_v_kmax > 1000);
  if (__pyx_2) {

    /* "/root/prime.pyx":10 */
    __pyx_v_kmax = 1000;
    goto __pyx_L2;
  }
  __pyx_L2:;

  /* "/root/prime.pyx":11 */
  __pyx_v_k = 0;

  /* "/root/prime.pyx":12 */
  __pyx_v_n = 2;

  /* "/root/prime.pyx":13 */
  while (1) {
    __pyx_L3:;
    __pyx_2 = (__pyx_v_k < __pyx_v_kmax);
    if (!__pyx_2) break;

    /* "/root/prime.pyx":14 */
    __pyx_v_i = 0;

    /* "/root/prime.pyx":15 */
    while (1) {
      __pyx_L5:;
      __pyx_2 = (__pyx_v_i < __pyx_v_k);
      if (__pyx_2) {
        __pyx_2 = ((__pyx_v_n % (__pyx_v_p[__pyx_v_i])) != 0);
      }
      if (!__pyx_2) break;

      /* "/root/prime.pyx":16 */
      __pyx_v_i = (__pyx_v_i + 1);
    }
    __pyx_L6:;

    /* "/root/prime.pyx":17 */
    __pyx_2 = (__pyx_v_i == __pyx_v_k);
    if (__pyx_2) {

      /* "/root/prime.pyx":18 */
      (__pyx_v_p[__pyx_v_k]) = __pyx_v_n;

      /* "/root/prime.pyx":19 */
      __pyx_v_k = (__pyx_v_k + 1);

      /* "/root/prime.pyx":20 */
      __pyx_1 = PyObject_GetAttr(__pyx_v_result, __pyx_n_append); if (!__pyx_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 20; goto __pyx_L1;}
      __pyx_3 = PyInt_FromLong(__pyx_v_n); if (!__pyx_3) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 20; goto __pyx_L1;}
      __pyx_4 = PyTuple_New(1); if (!__pyx_4) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 20; goto __pyx_L1;}
      PyTuple_SET_ITEM(__pyx_4, 0, __pyx_3);
      __pyx_3 = 0;
      __pyx_3 = PyObject_CallObject(__pyx_1, __pyx_4); if (!__pyx_3) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 20; goto __pyx_L1;}
      Py_DECREF(__pyx_1); __pyx_1 = 0;
      Py_DECREF(__pyx_4); __pyx_4 = 0;
      Py_DECREF(__pyx_3); __pyx_3 = 0;
      goto __pyx_L7;
    }
    __pyx_L7:;

    /* "/root/prime.pyx":21 */
    __pyx_v_n = (__pyx_v_n + 1);
  }
  __pyx_L4:;

  /* "/root/prime.pyx":22 */
  Py_INCREF(__pyx_v_result);
  __pyx_r = __pyx_v_result;
  goto __pyx_L0;

  __pyx_r = Py_None; Py_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1:;
  Py_XDECREF(__pyx_1);
  Py_XDECREF(__pyx_3);
  Py_XDECREF(__pyx_4);
  __Pyx_AddTraceback("prime.primes");
  __pyx_r = 0;
  __pyx_L0:;
  Py_DECREF(__pyx_v_result);
  return __pyx_r;
}

static __Pyx_InternTabEntry __pyx_intern_tab[] = {
  {&__pyx_n_append, "append"},
  {&__pyx_n_primes, "primes"},
  {0, 0}
};

static struct PyMethodDef __pyx_methods[] = {
  {"primes", (PyCFunction)__pyx_f_5prime_primes, METH_VARARGS|METH_KEYWORDS, 0},
  {0, 0, 0, 0}
};

static void __pyx_init_filenames(void); /*proto*/

PyMODINIT_FUNC initprime(void); /*proto*/
PyMODINIT_FUNC initprime(void) {
  __pyx_init_filenames();
  __pyx_m = Py_InitModule4("prime", __pyx_methods, 0, 0, PYTHON_API_VERSION);
  if (!__pyx_m) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5; goto __pyx_L1;};
  __pyx_b = PyImport_AddModule("__builtin__");
  if (!__pyx_b) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5; goto __pyx_L1;};
  if (PyObject_SetAttrString(__pyx_m, "__builtins__", __pyx_b) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5; goto __pyx_L1;};
  if (__Pyx_InternStrings(__pyx_intern_tab) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5; goto __pyx_L1;};

  /* "/root/prime.pyx":5 */
  return;
  __pyx_L1:;
  __Pyx_AddTraceback("prime");
}

static char *__pyx_filenames[] = {
  "prime.pyx",
};

/* Runtime support code */

static void __pyx_init_filenames(void) {
  __pyx_f = __pyx_filenames;
}

static int __Pyx_InternStrings(__Pyx_InternTabEntry *t) {
    while (t->p) {
        *t->p = PyString_InternFromString(t->s);
        if (!*t->p)
            return -1;
        ++t;
    }
    return 0;
}

#include "compile.h"
#include "frameobject.h"
#include "traceback.h"

static void __Pyx_AddTraceback(char *funcname) {
    PyObject *py_srcfile = 0;
    PyObject *py_funcname = 0;
    PyObject *py_globals = 0;
    PyObject *empty_tuple = 0;
    PyObject *empty_string = 0;
    PyCodeObject *py_code = 0;
    PyFrameObject *py_frame = 0;
    
    py_srcfile = PyString_FromString(__pyx_filename);
    if (!py_srcfile) goto bad;
    py_funcname = PyString_FromString(funcname);
    if (!py_funcname) goto bad;
    py_globals = PyModule_GetDict(__pyx_m);
    if (!py_globals) goto bad;
    empty_tuple = PyTuple_New(0);
    if (!empty_tuple) goto bad;
    empty_string = PyString_FromString("");
    if (!empty_string) goto bad;
    py_code = PyCode_New(
        0,            /*int argcount,*/
        0,            /*int nlocals,*/
        0,            /*int stacksize,*/
        0,            /*int flags,*/
        empty_string, /*PyObject *code,*/
        empty_tuple,  /*PyObject *consts,*/
        empty_tuple,  /*PyObject *names,*/
        empty_tuple,  /*PyObject *varnames,*/
        empty_tuple,  /*PyObject *freevars,*/
        empty_tuple,  /*PyObject *cellvars,*/
        py_srcfile,   /*PyObject *filename,*/
        py_funcname,  /*PyObject *name,*/
        __pyx_lineno,   /*int firstlineno,*/
        empty_string  /*PyObject *lnotab*/
    );
    if (!py_code) goto bad;
    py_frame = PyFrame_New(
        PyThreadState_Get(), /*PyThreadState *tstate,*/
        py_code,             /*PyCodeObject *code,*/
        py_globals,          /*PyObject *globals,*/
        0                    /*PyObject *locals*/
    );
    if (!py_frame) goto bad;
    py_frame->f_lineno = __pyx_lineno;
    PyTraceBack_Here(py_frame);
bad:
    Py_XDECREF(py_srcfile);
    Py_XDECREF(py_funcname);
    Py_XDECREF(empty_tuple);
    Py_XDECREF(empty_string);
    Py_XDECREF(py_code);
    Py_XDECREF(py_frame);
}