#!/usr/bin/env python

#This example script calculates r and kappa for FRET efficies from simulations with dyes.

import sys,os

from grompy import rvec
from grompy.meta import Trajectory
from grompy.index import rd_index

from numpy import array,savez,savetxt,float64
from numpy.linalg import norm

import warnings
warnings.simplefilter("ignore", RuntimeWarning)

def getPosVec(pos,ndx):
    a=array(pos[ndx[0]], dtype=float64)
    b=array(pos[ndx[1]], dtype=float64)
    c=array(pos[ndx[2]], dtype=float64)
    d=array(pos[ndx[3]], dtype=float64)
    pos=1./4*(a+b+c+d)
    vec=a+b-c-d
    vec=vec/norm(vec)
    return pos,vec

def getValues():
    isize,index,grpnames=rd_index(idxfile,2)
    if isize[0]!=4 or isize[1]!=4:
        print "Index sizes are not OK! Must be equal 4"
        sys.exit(1)

    traj=Trajectory(trajfile)
    values = []
    for positions,step,time,matrix,prec,box in traj:
        sys.stdout.write("Step: %d \tTime: %6.3f\r"%(step,time))
        donorpos,donorvec=getPosVec(positions,index[0])
        acceptorpos,acceptorvec=getPosVec(positions,index[1])
        R=norm(donorpos-acceptorpos)
        Rvec=(donorpos-acceptorpos)/R
        kappa=(donorvec*acceptorvec).sum()-3*((donorvec*Rvec).sum()*(acceptorvec*Rvec).sum())
        kappasq=kappa*kappa
        values.append((time,R,kappasq))

    return array(values)

print "Specify donor and acceptor index group with 4 atoms each."

if len(sys.argv)>1:
    trajfile=sys.argv[1]
else:
    trajfile="prot.xtc"
    
if len(sys.argv)>2:
    idxfile=sys.argv[2]
else:
    idxfile="donoracceptor.ndx"
    
if len(sys.argv)>3:
    outtextfile=sys.argv[3]
else:
    outtextfile="rkappa.dat"

if len(sys.argv)>4:
    outnpfile=sys.argv[4]
else:
    outnpfile="rkappa"

vals = getValues()
savez(outnpfile,vals)
savetxt(outtextfile,vals)
