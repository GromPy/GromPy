from ctypes import POINTER,c_int
from grompy.types import t_inputrec,XX,YY,ZZ

#from inputrec.h
def DEFORM(ir = t_inputrec):
    result = (ir.deform[XX][XX]!=0 or\
              ir.deform[YY][YY]!=0 or\
              ir.deform[ZZ][ZZ]!=0 or\
              ir.deform[YY][XX]!=0 or\
              ir.deform[ZZ][XX]!=0 or\
              ir.deform[ZZ][YY]!=0)
    return result
#define DEFORM(ir) ((ir).deform[XX][XX]!=0 || (ir).deform[YY][YY]!=0 || (ir).deform[ZZ][ZZ]!=0 || (ir).deform[YY][XX]!=0 || (ir).deform[ZZ][XX]!=0 || (ir).deform[ZZ][YY]!=0)