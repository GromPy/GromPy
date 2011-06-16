from ctypes import c_char_p,c_int,pointer,POINTER,byref
from grompy import libgmx


def rd_index(statfile,ngroups,isizep=None,indexp=None,grpnamesp=None):
    filename=c_char_p(statfile)
    if not isizep:
        isize=(c_int*ngroups)()
        isizep=pointer(isize)
    if not indexp:
        index=(POINTER(c_int)*ngroups)()
        indexp=pointer(index)
    if not grpnamesp:  
        grpnames=(c_char_p*ngroups)()
        grpnamesp=pointer(grpnames)
    
    libgmx.rd_index(filename,ngroups,isizep,indexp,grpnamesp)
  
    return isize,index,grpnames
  
  
def get_index(atoms,ngroups,filename=c_char_p(),isizep=None,indexp=None,grpnamesp=None):
    if type(filename)==str:
        filename=c_char_p(filename)

    if not isizep:
        isize=(c_int*ngroups)()
        isizep=pointer(isize)
    if not indexp:
        index=(POINTER(c_int)*ngroups)()
        indexp=pointer(index)
    if not grpnamesp:  
        grpnames=(c_char_p*ngroups)()
        grpnamesp=pointer(grpnames)

    libgmx.get_index(byref(atoms),filename,ngroups,isizep,indexp,grpnamesp)
    return isize,index,grpnames

