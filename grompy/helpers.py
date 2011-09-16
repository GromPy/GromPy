def rvecsub(a,b):
    for i in range(3):
        a[i]=a[i]-b[i]

def rvecadd(a,b):
    for i in range(3):
        a[i]=a[i]+b[i]

def rvecscale(a,b):
    for i in range(3):
        a[i]=a[i]*b
                
def rvecset(a,b):
    for i in range(3):
        a[i]=b[i]
        
def coords2PDB(coords,pdbfname,AtomName="X",AtomAltLoc="",AtomResName="",AtomChainID="",AtomResSeq=1,AtomICode="",AtomOccupancy=0.0,AtomTempFactor=0.0,AtomSegID="",AtomElement="",AtomCharge="",AtomSuffix=""):
    fh=open(pdbfname,"w")
    for i in range(len(coords)):
        fh.write("ATOM  %5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s%s\n" % (i+1,AtomName,AtomAltLoc,AtomResName,AtomChainID,AtomResSeq,AtomICode,coords[i][0],coords[i][1],coords[i][2],AtomOccupancy,AtomTempFactor,AtomSegID,AtomElement,AtomCharge,AtomSuffix))
    fh.close()
    print pdbfname,"written."

