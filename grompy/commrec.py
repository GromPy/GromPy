from ctypes import POINTER,c_int
from grompy.types import t_commrec, gmx_multisim_t

def MASTERNODE(cr = POINTER(t_commrec)):
    return (c_int(cr.contents.nodeid).value == c_int(0).value)

def MASTERTHREAD(cr = POINTER(t_commrec)):
    return (c_int(cr.contents.threadid).value == c_int(0).value)

def MASTER(cr = POINTER(t_commrec)):
    return (MASTERNODE(cr) and MASTERTHREAD(cr))
#from commrec.h
##define MASTERNODE(cr)     ((cr)->nodeid == 0)
##define MASTERTHREAD(cr)   ((cr)->threadid == 0)
##define MASTER(cr)         (MASTERNODE(cr) && MASTERTHREAD(cr))

#from commrec.h
def NODEPAR(cr = POINTER(t_commrec)):
    return (c_int(cr.contents.nnodes).value > c_int(1).value)
#define NODEPAR(cr)        ((cr)->nnodes > 1)
def THREADPAR(cr = POINTER(t_commrec)):
    return (c_int(cr.contents.nthreads).value > c_int(1).value)
#define THREADPAR(cr)      ((cr)->nthreads > 1)
def PAR(cr = POINTER(t_commrec)):
    return (NODEPAR(cr) or THREADPAR(cr))
#define PAR(cr)            (NODEPAR(cr) || THREADPAR(cr))

#from commrec.h
def MULTISIM(cr = POINTER(t_commrec)):
    return bool(cr.contents.ms)
#define MULTISIM(cr)       ((cr)->ms)
def MASTERSIM(ms = POINTER(gmx_multisim_t)):
    if(bool(ms)):
        return (not bool(ms))
    else:
        return (not bool(ms.contents.sim))
#define MASTERSIM(ms)      ((ms)->sim == 0)
def MULTIMASTER(cr = POINTER(t_commrec)):
    return (MASTER(cr) and (not MULTISIM(cr) or MASTERSIM(cr.contents.ms)))
#/* The master of all (the node that prints the remaining run time etc.) */
#define MULTIMASTER(cr)    (MASTER(cr) && (!MULTISIM(cr) || MASTERSIM((cr)->ms)))

#from commrec.h:
DUTY_PP  = (1<<0)
DUTY_PME = (1<<1)
#define DUTY_PP  (1<<0)
#define DUTY_PME (1<<1)
def SIMMASTER(cr = POINTER(t_commrec)):
    return (MASTER(cr) and (cr.contents.duty & DUTY_PP))
#define SIMMASTER(cr)      (MASTER(cr) && ((cr)->duty & DUTY_PP))

def DOMAINDECOMP(cr = POINTER(t_commrec)):
    return (bool(cr.contents.dd))
#define DOMAINDECOMP(cr)   ((cr)->dd != NULL)