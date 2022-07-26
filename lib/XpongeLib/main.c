#include "parmchk_mod/parmchk2.c"

#include <Python.h>

PyObject* _parmchk2(PyObject* self, PyObject* args)
{
    char *i;
    char *iformat;
    char *o;
    char *datapath;
    int print_all, print_dihedral_contain_X, gaffORgaff2;
    if (!PyArg_ParseTuple(args, "ssssiii", &i, &iformat, &o, &datapath, &print_all, &print_dihedral_contain_X, &gaffORgaff2))
    {
        return NULL;
    }
    parmchk2(i, iformat, o, datapath, print_all, print_dihedral_contain_X, gaffORgaff2);
    return Py_None;
}

static PyMethodDef XpongeLibMethods[] = {
{"_parmchk2", _parmchk2, METH_VARARGS, ""},
{NULL, NULL, 0, NULL}
};

static struct PyModuleDef XpongeLibModule = {
PyModuleDef_HEAD_INIT,
"XpongeLib",
"the C++ Lib for the package for building molecular dynamics inputs for SPONGE",
-1,
XpongeLibMethods
};

PyMODINIT_FUNC PyInit_backend(void)
{
    return PyModule_Create(&XpongeLibModule);
}
