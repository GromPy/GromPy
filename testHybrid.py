#from ctypes import c_char_p
import sys

TprDir  = sys.argv[1]
Mu      = float(sys.argv[2])
NStart  = int(sys.argv[3])
#NMax    = 426
#NMin    = 0
MolName = sys.argv[4]
NMin    = int(sys.argv[5])
NMax    = int(sys.argv[6])
if(len(sys.argv)>7 and sys.argv[7]=='NrgGrps'):
    print "This implementation does not work!!"
    print "It needs some extra coding: for each GCMC move an energy"
    print "evaluation has to be performed IF the previous GCMC move"
    print "was rejected."
    print "Such a scheme will result in the same amount of energy"
    print "evaluations as the other implementation."
    print "EXITING..."
    sys.exit(1)
    bUseNrgGrps = True
    from grompy.HybridNrgGrps import *
    hybrid_nrggrps(Mu,NStart,NMax,NMin,TprDir)
else:
    bUseNrgGrps = False
    from grompy.Hybrid import *
    hybrid(Mu,NStart,NMax,NMin,TprDir,MolName)
