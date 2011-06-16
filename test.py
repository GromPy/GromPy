from grompy.fio import *
from grompy.tpxio import *
from grompy.index import *

print "opening trajectory file"
fh = open_xtc("test/traj.xtc","r")

print "reading first frame"
natoms,step,time,matrix,vecpnt,prec,bok,ret = read_first_xtc(fh)
print "natoms",natoms,"step",step,"time",time,"box",matrix,"vecpnt",vecpnt,"prec",prec,"bok",bok

print "first atoms"
for i in xrange(10):
	print vecpnt[i][0],vecpnt[i][1],vecpnt[i][2]

step,time,matrix,prec,bok,ret = read_next_xtc(fh,natoms,vecpnt)
print "step",step,"time",time,"box",matrix,"prec",prec,"bok",bok

print "first atoms"
for i in xrange(10):
        print vecpnt[i][0],vecpnt[i][1],vecpnt[i][2]

print "closing file"
close_xtc(fh)

print "opening topology"

title,top,ePBC,xp,vp,box,bMass = read_tps_conf("test/topol.tpr")
print "title",title,"topology",top,"ePBC",ePBC,"xpointer",xp,"vpointer",vp,"box",box,"bmass",bMass

print "ok now reading indices"
isize1,index1,grpnames1 = rd_index("test/index.ndx",2)
print "without index file"
isize2,index2,grpnames2 = get_index(top.atoms,1)
print "with index file"
isize3,index3,grpnames3 = get_index(top.atoms,1,filename="test/index.ndx")
print "ok - now printing the indices"

print isize1[0],index1[0][0],grpnames1[0]
print isize2[0],index2[0][0],grpnames2[0]
print isize3[0],index3[0][0],grpnames3[0]
print "Done"
