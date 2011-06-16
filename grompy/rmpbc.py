from grompy import libgmx

def rm_pbc(idefpointer,natoms,box,xpointer,ePBC=0,xoutpointer=None):
    if not xoutpointer:
        xoutpointer=xpointer
    libgmx.rm_pbc(idefpointer,ePBC,natoms,box,xpointer,xoutpointer)
    
def rm_gropbc(atoms,xpointer,box):
    libgmx.rm_gropbc(atoms,xpointer,box)
