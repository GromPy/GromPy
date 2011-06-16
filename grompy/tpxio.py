from ctypes import c_char,byref,POINTER,c_int
from grompy.types import t_topology
from grompy import libgmx,GMXctypesError,rvec,matrix

def read_tps_conf(filename):
  f = tpxfile(filename)
  return f.title,f.top,f.ePBC,f.x,f.v,f.box,f.bMass
  
class tpxfile:
    def __init__(self,filename):
        toptitle=(c_char*256)()
        top=t_topology()
        xp=POINTER(rvec)()
        vp=POINTER(rvec)()
        cePBC=c_int(-1)
        box=matrix()
        bMass=c_int()
        ret = libgmx.read_tps_conf(filename,byref(toptitle),byref(top),byref(cePBC),byref(xp),byref(vp),box,byref(bMass))
        if ret!=1:
            raise GMXctypesError,"Error reading topology with read_tps_conf"
        self.title,self.top,self.ePBC,self.x,self.v,self.box,self.bMass=toptitle.value,top,cePBC.value,xp,vp,box,bMass.value
    def __iter__(self):
        for i in range(self.top.atoms.nr):
            yield atom(self,i)
    def __getitem__(self,i):
            return atom(self,i)
          
class atom:
  def __init__(self,tpx,nr):
    self.tpx,self.nr=tpx,nr
    self.atom=self.tpx.top.atoms.atom[self.nr]
  
  #wrapper can be written with either @property syntax or __getattr__
  #using getattr for those attributes where several are similar
  def __getattr__(self,name):
    if name in ["m","q","chain","resnr","type"]:
      return getattr(self.atom,name)
    elif name in ["x","v"]:
      return getattr(self.tpx,name)[self.nr]
    elif name in ["atomname","atomtype"]:
      return getattr(self.tpx.top.atoms,name)[self.nr][0]
    elif name in ["radius","vol","surftens","atomnumber"]:
      return getattr(self.tpx.top.atomtypes,name)[self.type]
    else:
      raise AttributeError(name)
  #using property syntax for odd ones
  @property
  def resname(self): return self.tpx.top.atoms.resname[self.resnr][0]
