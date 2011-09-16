from grompy.fio import read_first_xtc,open_xtc,close_xtc,read_next_xtc

class Frame:
    def __init__(self,vecpnt,step,time,matrix,prec,bok):
        self.coords=vecpnt
        self.step=step
        self.time=time
        self.box=matrix
        self.precision=prec
        self.bok=bok

class Trajectory:
    def __init__(self,trajname,retframes=False):
        self.trajname=trajname
        print "Opening Trajectory %s"%self.trajname
        self.fh = open_xtc(trajname,"r")
        self.natoms,self.step,self.time,self.matrix,self.vecpnt,self.prec,self.bok,retcode = read_first_xtc(self.fh)
        self.first=True
        self.retframes=retframes
        if not retcode==1:
            raise StopIteration
        
    def __iter__(self):
        return self
    
    
    def next(self):
        if self.first:
            self.first=False
            if not self.retframes:
                return self.vecpnt,self.step,self.time,self.matrix,self.prec,self.bok
            return Frame(self.vecpnt,self.step,self.time,self.matrix,self.prec,self.bok)
        
        else:
            self.step,self.time,self.matrix,self.prec,self.bok,retcode = read_next_xtc(self.fh,self.natoms,self.vecpnt)
            if not retcode==1:
                raise StopIteration
            if not self.retframes:
                return self.vecpnt,self.step,self.time,self.matrix,self.prec,self.bok
            return Frame(self.vecpnt,self.step,self.time,self.matrix,self.prec,self.bok)
        
    def __del__(self):
        print "Closing trajectory %s"%self.trajname
        close_xtc(self.fh)
        


def rvec2str(rvec):
    return "(%f, %f, %f)" % (rvec[0],rvec[1],rvec[2])

def rvecs2str(rvecs,count=0):
    if count==0:
        count=len(rvecs)
    strg=""
    for i in range(count):
        strg+="%d:"%(i+1)
        strg+=rvec2str(rvecs[i])
        strg+="\n"
    return strg