from ctypes import c_char_p,POINTER,c_int,byref
from grompy import libgmx,rvec,matrix
from grompy import c_real

def open_xtc(filename,mode):
	fileptr=c_char_p(filename)
	modeptr=c_char_p(mode)
	return libgmx.open_xtc(fileptr,modeptr)
	
def close_xtc(handle):
	libgmx.close_xtc(handle)

def read_first_xtc(handle):
	natoms = c_int()
	step = c_int()
	time = c_real()
	prec = c_real()
	bOK = c_int()
	box =  matrix()
	xpointer = POINTER(rvec)()
	ret=libgmx.read_first_xtc(handle,byref(natoms),byref(step),byref(time),box,byref(xpointer),byref(prec),byref(bOK))
#	if ret!=1:
#		raise GMXctypesError,"read_first_xtc did not return 1"
	return natoms.value,step.value,time.value,box,xpointer,prec.value,bOK.value,ret
	
def read_next_xtc(handle, natoms, xpointer):
	step = c_int()
	time = c_real()
	prec = c_real()
	bOK = c_int()
	box =  matrix()
	ret=libgmx.read_next_xtc(handle,natoms,byref(step),byref(time),box,xpointer,byref(prec),byref(bOK))
#	if ret!=1:
#		raise GMXctypesError,"read_next_xtc did not return 1"

	return step.value,time.value,box,prec.value,bOK.value,ret

def write_xtc(handle,natoms,step,time,box,xpointer,prec):
	ret=libgmx.write_xtc(handle,natoms,step,time,box,xpointer,prec)
	return ret
#	if ret!=1:
#		raise GMXctypesError,"write_xtc did not return 1"
	
