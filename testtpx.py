from grompy.fio import *
from grompy.tpxio import *
from grompy.index import *
from itertools import *


print "opening topology"

#testing backward compatible function
title,top,ePBC,xp,vp,box,bMass = read_tps_conf("test/topol.tpr")
print "title",title,"topology",top,"ePBC",ePBC,"xpointer",xp,"vpointer",vp,"box",box,"bmass",bMass

#new class syntax
f = tpxfile("test/topol.tpr")
#title, PBC, box, bMass
print "title",f.title,"ePBC",f.ePBC,"box",f.box,"bmass",f.bMass,
#topology, coords, vel: can still be accessed directly - better atom interface - see next
print "topology",f.top,"xpointer",f.x,"vpointer",f.v

#iterate over atoms
#for a in islice(f,10):  #first 10, all without islice
for a in f:
  print a.x[0],a.x[1],a.x[2],a.m,a.atomname,a.atomtype,a.type,a.resname,a.q,a.chain,a.radius
#you can also access atoms by index
print f[9].m   #mass just as example


