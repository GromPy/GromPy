from ctypes import c_int,POINTER
from grompy import libgmx,matrix
from grompy import c_real

libgmx.calc_similar_ind.restype = c_real


def rmsdev(natoms,indices,masses,atomsvec1,atomsvec2):
    return libgmx.calc_similar_ind(c_int(0),natoms,indices,masses,atomsvec1,atomsvec2)

def rhodev(natoms,indices,masses,atomsvec1,atomsvec2):
    return libgmx.calc_similar_ind(c_int(1),natoms,indices,masses,atomsvec1,atomsvec2)
    
def calc_fit_R(natoms,weights,atomsvec1,atomsvec2,dimensions=3):
    mtx = matrix()
    libgmx.calc_fit_R(dimensions,natoms,weights,atomsvec1,atomsvec2,mtx)
    return mtx

def do_fit(natoms,weights,atomsvec1,atomsvec2,dimensions=3):
    libgmx.do_fit_ndim(dimensions,natoms,weights,atomsvec1,atomsvec2)
    
def reset_x(ncmatoms,cmindices,atomsvec,masses,nrtoreset,dimensions=3,rindices=None):
    if not rindices:
        rindices=POINTER(c_int)()
    libgmx.reset_x_ndim(dimensions,ncmatoms,cmindices,nrtoreset,rindices,atomsvec,masses)
    