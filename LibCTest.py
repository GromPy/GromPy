import sys,unittest
from ctypes import cdll,c_char_p,CFUNCTYPE,c_int,Structure,addressof
from ctypeslib.contrib.pythonhdr import *

libc = cdll.LoadLibrary("/lib/libc.so.6")

s = c_char_p("G R O M A C S")
#s = "G R O M A C S"
libc.printf("%s\n",s)

class POINT(Structure):
	_fields_ = [("x", c_int), ("y", c_int)]

TenPointsArrayType = POINT * 2 
arr = TenPointsArrayType(POINT(x=c_int(0),y=1),
                         POINT(x=2,y=3))
for pt in arr:
    print pt.x, pt.y, addressof(pt)

def test_file(self):
        closefn = CFUNCTYPE(c_int, FILE_ptr)
        def close(f):
            return 1
        close_callback = closefn(close)
        try:
            path = sys.executable
            f = file(path, 'rb')
            try:
                fp = PyFile_AsFile(f)
                self.failUnless(fp)
                g = PyFile_FromFile(fp, path, 'rb', close_callback)
                fno = g.fileno()
                del g
                self.failUnlessEqual(f.fileno(), fno)
            finally:
                f.close()
        except (NameError, IOError):
            pass

print test_file(joehoe)
