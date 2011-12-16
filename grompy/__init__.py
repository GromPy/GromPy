from ctypes import c_float,\
                   cdll,\
                   c_double,\
                   c_int,\
                   pythonapi,\
                   py_object,\
                   c_void_p

from os import environ

import sys

libcname="/lib/libc.so.6" # default path, but is different in e.g. Ubuntu 11.10
if environ.has_key("LIBCPATH"):
    libcname=environ["LIBCPATH"]

if environ.has_key("GROMPY_DOUBLE"):
    isdouble=True
    c_real = c_double
    print "Loading GromPy with double precision library..."
    libmdname="libmd_d.so"
    libgmxname="libgmx_d.so"
    libmdrunname="mdrun_d.so"

else:
    isdouble=False
    c_real   = c_float
    print "Loading GromPy with single precision library..."
    libmdname="libmd.so"
    libgmxname="libgmx.so"
    libmdrunname="mdrun.so"
    
if environ.has_key("GROMPY_LIBEXT"):
    libmdname="libmd%s.so"%environ["GROMPY_LIBEXT"]
    libgmxname="libgmx%s.so"%environ["GROMPY_LIBEXT"]
    libmdrunname="mdrun%s.so"%environ["GROMPY_LIBEXT"]
    
if environ.has_key("GROMPY_LIBGMX"):
    libgmxname=environ["GROMPY_LIBGMX"]

if environ.has_key("GROMPY_LIBMD"):
    libmdname=environ["GROMPY_LIBMD"]

if environ.has_key("GROMPY_LIBMDRUN"):
    libmdrunname=environ["GROMPY_LIBMDRUN"]

libc     = cdll.LoadLibrary(libcname)
libmd    = cdll.LoadLibrary(libmdname)
libgmx   = cdll.LoadLibrary(libgmxname)
libmdrun = cdll.LoadLibrary(libmdrunname)

# gromacs shuold be compiled with -lfftw3f!!!
#    libmdrun = cdll.LoadLibrary(environ["HOME"]+"/src/gromacs-4.0.5_TEST/src/kernel/libmdrun0/mdrun.so") # for debugging purposes


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
