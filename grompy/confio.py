from ctypes import c_char_p,byref,POINTER,c_int
from grompy.types import atom_id
from grompy import libgmx,matrix

def init_t_atoms(atoms,natoms,bPDBinfo=False):
    if bPDBinfo==True:
        pdbinf = c_int(1)
    else:
        pdbinf = c_int(0)
        
    libgmx.init_t_atoms(byref(atoms),natoms,pdbinf)

def read_g96_conf(fhandle,filename,trxframe):
    fname=c_char_p(filename)
    libgmx.read_g96_conf(fhandle,fname,byref(trxframe))
    

def write_g96_conf(fhandle,trxframe,nindex,index=None):
    if not index:
        ndx=POINTER(atom_id)()
    libgmx.write_g96_conf(fhandle,byref(trxframe),nindex,ndx)

def gro_next_x_or_v(fhandle,trxframe):
    libgmx.gro_next_x_or_v(fhandle,byref(trxframe))

def gro_first_x_or_v(fhandle,trxframe):
    return libgmx.gro_first_x_or_v(fhandle,byref(trxframe))

def write_hconf_p(fhandle,title,atoms,ndec,x,v,box,nindex=None,index=None):
    ttl=c_char_p(title)
    if nindex:
        libgmx.write_hconf_indexed_p(fhandle,ttl,byref(atoms),nindex,index,ndec,x,v,box)
    else:
        libgmx.write_hconf_p(fhandle,ttl,byref(atoms),ndec,x,v,box)

def write_sto_conf(outfile,title,atoms,x,v,ePBC,box,nindex=None,index=None):
    ttl=c_char_p(title)
    outf=c_char_p(outfile)
    print outf,ttl,byref(atoms),x,v,ePBC,box,nindex,index
    if nindex:
        print "writing"
        libgmx.write_sto_conf_indexed(outf,ttl,byref(atoms),x,v,ePBC,box,nindex,index)
    else:
        libgmx.write_sto_conf(outf,ttl,byref(atoms),x,v,ePBC,box)

def write_sto_conf_mtop(outfile,title,mtop,x,v,ePBC,box):
    ttl=c_char_p(title)
    outf=c_char_p(outfile)
    libgmx.write_sto_conf_mtop(outf,ttl,byref(mtop),x,v,ePBC,box)

def get_stx_coordnum(infile):
    natoms=c_int(-1)
    ifile=c_char_p(infile)
    libgmx.get_stx_coordnum(ifile,byref(natoms))
    return natoms

def read_stx_conf(infile,title,atoms,x,v):
    ifile=c_char_p(infile)
    ttl=c_char_p(title)
    m=matrix()
    pbc=c_int(-1)
    libgmx.read_stx_conf(ifile, ttl,byref(atoms),x,v,byref(pbc),m);
    return pbc,m
