from ctypes import c_char_p

from grompy.md import *

filename = "test.log"
fplog    = c_char_p(filename)

mdrun()