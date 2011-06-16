from ctypes import c_double
def difftime(end,start):
    return c_double(float(end.value)-float(start.value))
#from sim_util.c
#define difftime(end,start) ((double)(end)-(double)(start))