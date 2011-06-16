from ctypes import c_float,\
                   cdll,\
                   c_double,\
                   c_int,\
                   pythonapi,\
                   py_object,\
                   c_void_p

from os import environ

import sys

libc = cdll.LoadLibrary("/lib/libc.so.6")
# gromacs shuold be compiled with -lfftw3f!!!
if environ.has_key("GROMPYDOUBLE"):
    isdouble=True
    print "Loading GromPy with double precision library..."
    c_real = c_double
    libmd  = cdll.LoadLibrary("libmd_d.so")
    libgmx = cdll.LoadLibrary("libgmx_d.so")
else:
    isdouble=False
    print "Loading GromPy with single precision library..."
    c_real   = c_float
    libmd    = cdll.LoadLibrary("libmd.so")
    libgmx   = cdll.LoadLibrary("libgmx.so")
#    libmdrun = cdll.LoadLibrary(environ["HOME"]+"/src/gromacs-4.0.5_TEST/src/kernel/libmdrun0/mdrun.so") # for debugging purposes
    libmdrun = cdll.LoadLibrary(environ["HOME"]+"/src/gromacs-4.0.5_TEST/src/kernel/libmdrun1/mdrun.so")

#FILE * to stderr, stdout
pythonapi.PyFile_AsFile.restype = c_void_p
stderr    = pythonapi.PyFile_AsFile(py_object(sys.stderr))
stdout    = pythonapi.PyFile_AsFile(py_object(sys.stdout))

rvec      = c_real*3
dvec      = c_double*3
ivec      = c_int*3
splinevec = c_real*3
matrix    = c_real*3*3
tensor    = c_real*3*3

class GMXctypesError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
