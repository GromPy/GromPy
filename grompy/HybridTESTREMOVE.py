import sys
import os
import math
import random
import re
from scipy import array,rand,dot
from ctypes import c_char_p,\
                   POINTER,\
                   c_int,\
                   addressof,\
                   byref,\
                   c_ulong,\
                   c_char
from grompy import libmdrun,rvec,c_real,ivec,stderr,libc,stdout,tensor
import grompy.types as gt
import grompy.commrec as commrec
import grompy.statutil as statutil
import grompy.physics as phys
import grompy.RandomRotation as randrot
from types import mdrunner_comm

iEps    = 1.0/5.0
TStar   = 5.0*phys.KILO/phys.RGAS
TRed    = 900.0/TStar
BetaRed = 1.0/TRed

def GetAtomTypesOfMolecule(mc = mdrunner_comm, MolType = str):
    AtomType = []

    iMolTypeIndex = 0
    for i in range(mc.mtop.contents.nmoltype):
        if(mc.mtop.contents.moltype[i].name.contents.value==MolType):
            iMolTypeIndex = i
            break
    NAtomsPerMolType = mc.mtop.contents.moltype[iMolTypeIndex].atoms.nr

    for i in range(NAtomsPerMolType):
        AtomType.append(mc.mtop.contents.moltype[iMolTypeIndex].atoms.atom[i].type)

    return AtomType

def GetIndexOffSet(mc = mdrunner_comm, MolType = str):
    NAtomsPerMolType = []
    for i in range(mc.mtop.contents.nmoltype):
        NAtomsPerMolType.append(mc.mtop.contents.moltype[i].atoms.nr)

    NMolPerMolType   = []
    for i in range(mc.mtop.contents.nmolblock):
        NMolPerMolType.append(mc.mtop.contents.molblock[i].nmol)

    MolTypeIndex = 0
    for i in range(mc.mtop.contents.nmoltype):
        if(mc.mtop.contents.moltype[i].name.contents.value==MolType):
            MolTypeIndex = i
            break

    MolBlockIndex = 0
    for i in range(mc.mtop.contents.nmolblock):
        if(mc.mtop.contents.molblock[i].type==MolTypeIndex):
            MolBlockIndex = i
            break

    IndexOffSet = 0
    for i in range(MolBlockIndex):
        IndexOffSet += i*NMolPerMolType[i]*NAtomsPerMolType[i]

    return IndexOffSet, MolBlockIndex

def RotateMoleculeRandomly(Conf=[]):
#    rvecarray = rvec*NAtomsOfMolType
    ConfNew   = []
    randrotm  = randrot.RandRotation(rand(3))
    for j in range(len(Conf)):
        x     = array([Conf[j][gt.XX],Conf[j][gt.YY],Conf[j][gt.ZZ]])
        randx = dot(randrotm,x)
        ConfNew.append([randx[gt.XX],randx[gt.YY],randx[gt.ZZ]])
#        ConfNew[j][gt.XX] = c_real(randx[gt.XX])
#        ConfNew[j][gt.YY] = c_real(randx[gt.YY])
#        ConfNew[j][gt.ZZ] = c_real(randx[gt.ZZ])
        del x, randx
    return ConfNew

def GetFirstMoleculeAndPutFirstAtomAtOrigin(mc=gt.mdrunner_comm,MolType=str):
    NAtomsPerMolType = []
    for i in range(mc.mtop.contents.nmoltype):
        NAtomsPerMolType.append(mc.mtop.contents.moltype[i].atoms.nr)

    NMolPerMolType   = []
    for i in range(mc.mtop.contents.nmolblock):
        NMolPerMolType.append(mc.mtop.contents.molblock[i].nmol)

    MolTypeIndex = 0
    for i in range(mc.mtop.contents.nmoltype):
        if(mc.mtop.contents.moltype[i].name.contents.value==MolType):
            MolTypeIndex = i
            break

    MolBlockIndex = 0
    for i in range(mc.mtop.contents.nmolblock):
        if(mc.mtop.contents.molblock[i].type==MolTypeIndex):
            MolBlockIndex = i
            break

    NAtomsOfMolType = mc.mtop.contents.moltype[MolTypeIndex].atoms.nr
    IndexOffSet = 0
    for i in range(MolBlockIndex):
        IndexOffSet += i*NMolPerMolType[i]*NAtomsPerMolType[i]

#    rvecarray = rvec*NAtomsOfMolType
    Conf      = []
    for i in range(NAtomsOfMolType):
        Conf.append([(mc.state.contents.x[i+IndexOffSet][gt.XX]-mc.state.contents.x[IndexOffSet+0][gt.XX]),
                     (mc.state.contents.x[i+IndexOffSet][gt.YY]-mc.state.contents.x[IndexOffSet+0][gt.YY]),
                     (mc.state.contents.x[i+IndexOffSet][gt.ZZ]-mc.state.contents.x[IndexOffSet+0][gt.ZZ])])
    return Conf

def CheckIfTrialMolIsFirstInTopology(mc=gt.mdrunner_comm,iTrialMol=int,MolBlock=int):
    gtypes = (c_char_p*gt.egcNR.value).in_dll(libmdrun,"gtypes")
    gtypes_ener_index = 0
    for i in range(gt.egcNR.value):
        if(gtypes[i]=="Energy Mon."):
            gtypes_ener_index = i
            break

    NAtomsMol = mc.mtop.contents.molblock[MolBlock].natoms_mol
    iTrialRemoveAtomList = []
    for i in range(NAtomsMol):
        if(bool(mc.mtop.contents.groups.grpnr[gtypes_ener_index])):
            if(mc.mtop.contents.groups.grpnr[gtypes_ener_index][i]==0):
                iTrialRemoveAtomList.append(i)

    bFirstMoleculeInList = False
    for i in iTrialRemoveAtomList:
        if(i==0):
            bFirstMoleculeInList = True
    if(not bFirstMoleculeInList):
        print "Trial removal particle is NOT the first in the list. Regenerate the index file(s)!!!"
        sys.exit(1)

    return

def SelectMoleculeForTrialRemoval(mc=gt.mdrunner_comm,CurrNMol=int,MolBlock=int,iTrialMol=int):
    iMolNewRemove = random.randint(0,CurrNMol-1)
    libmdrun.swap_coordinate_indices(byref(mc),
                                     c_int(MolBlock),
                                     c_int(iTrialMol),
                                     c_int(iMolNewRemove))
    return iMolNewRemove

def WriteCoords(file=str,mc=gt.mdrunner_comm,n=int):
    fw=open(file,'w')
    fw.write(str(n)+'\n')
    fw.write('coordinates\n')
    for i in range(n):
        x  = mc.state.contents.x[i]
        fw.write('O '+str(10.0*x[gt.XX])+' '+str(10.0*x[gt.YY])+' '+str(10.0*x[gt.ZZ])+'\n')
    fw.close()
    vlj  = 0.0
    s2   = 0.47*0.47
    box  = mc.box[gt.XX][gt.XX]
    hbx  = box/2.0
    rc2  = 1.2*1.2
    ecut = 4.0*5.0*(math.pow((s2/rc2),6.0)-math.pow((s2/rc2),3.0))
    for i in range(0,n-1):
        xi  = mc.state.contents.x[i]
        for j in range(i+1,n):
            if(i!=j):
                xj  = mc.state.contents.x[j]
                dx  = xj[gt.XX] - xi[gt.XX]
                if(dx>hbx):
                    dx -= box
                elif(dx<-hbx):
                    dx += box
                dy  = xj[gt.YY] - xi[gt.YY]
                if(dy>hbx):
                    dy -= box
                elif(dy<-hbx):
                    dy += box
                dz  = xj[gt.ZZ] - xi[gt.ZZ]
                if(dz>hbx):
                    dz -= box
                elif(dx<-hbx):
                    dz += box
                dr2 = dx*dx + dy*dy + dz*dz
                if(dr2<rc2):
                    vlj += 5.0*4.0*(math.pow((s2/dr2),6.0)-math.pow((s2/dr2),3.0))-ecut
    print vlj
    return

def GetDeBroglieLambdaIn_nm(Mass = float,Temperature = float):
    Lambda = math.sqrt(math.pow(phys.PLANCK1,2.0)/ \
                       (2.0*phys.M_PI*Mass*phys.AMU*phys.BOLTZMANN*Temperature))/ \
             phys.NANO
    return Lambda

def GetTemperatureOfTpr(mc = gt.mdrunner_comm):
    return mc.inputrec.contents.opts.ref_t.contents.value

def MCInsert(currmc = gt.mdrunner_comm,
             trialmc = gt.mdrunner_comm,
             Ncurr = int,
             Ntrial = int,
             NAtomsOfMolType = int,
             SigmaVel = [],
             SigmasAtomCrossTerms2Frac = [[]],
             IndexOffSet = int,
             AtomType = [],
             Conf = []):

    c_intarray = c_int*NAtomsOfMolType
    c_atomtype = c_intarray()
    for i in range(NAtomsOfMolType):
        c_atomtype[i] = c_int(AtomType[i])

    c_realarray          = c_real*(len(SigmasAtomCrossTerms2Frac)*len(SigmasAtomCrossTerms2Frac))
    SigmaCrossTerms2Frac = c_realarray()
    for i in range(len(SigmasAtomCrossTerms2Frac)):
        for j in range(len(SigmasAtomCrossTerms2Frac)):
            index = i*len(SigmasAtomCrossTerms2Frac)+j
            SigmaCrossTerms2Frac[index] = c_real(SigmasAtomCrossTerms2Frac[i][j])

    rvecarray = rvec*NAtomsOfMolType
    vecx      = rvecarray()
    vecv      = rvecarray()
    r         = [0.0]*3
    r[gt.XX]  = random.random()*trialmc.box[gt.XX][gt.XX]
    r[gt.YY]  = random.random()*trialmc.box[gt.YY][gt.YY]
    r[gt.ZZ]  = random.random()*trialmc.box[gt.ZZ][gt.ZZ]
    ConfNew   = []
    if(NAtomsOfMolType>1):
        ConfNew = RotateMoleculeRandomly(Conf)
    else:
        ConfNew = Conf
#   maybe COM vel + random rotations?
    for i in range(NAtomsOfMolType):
        vecx[i][gt.XX] = r[gt.XX] + ConfNew[i][gt.XX]
        vecx[i][gt.YY] = r[gt.YY] + ConfNew[i][gt.YY]
        vecx[i][gt.ZZ] = r[gt.ZZ] + ConfNew[i][gt.ZZ]
        vecv[i][gt.XX] = random.gauss(0.0,SigmaVel[i])
        vecv[i][gt.YY] = random.gauss(0.0,SigmaVel[i])
        vecv[i][gt.ZZ] = random.gauss(0.0,SigmaVel[i])

    xadd = rvecarray()
    vadd = rvecarray()
    for i in range(NAtomsOfMolType):
        xadd[i][gt.XX] = c_real(vecx[i][gt.XX])
        xadd[i][gt.YY] = c_real(vecx[i][gt.YY])
        xadd[i][gt.ZZ] = c_real(vecx[i][gt.ZZ])
        vadd[i][gt.XX] = c_real(vecv[i][gt.XX])
        vadd[i][gt.YY] = c_real(vecv[i][gt.YY])
        vadd[i][gt.ZZ] = c_real(vecv[i][gt.ZZ])

    trialmc.inputrec.contents.nsteps    = c_int(0)
    trialmc.inputrec.contents.nstlog    = c_int(0)
    trialmc.inputrec.contents.nstxout   = c_int(0)
    trialmc.inputrec.contents.nstvout   = c_int(0)
    trialmc.inputrec.contents.nstfout   = c_int(0)
    trialmc.inputrec.contents.nstenergy = c_int(0)
    trialmc.inputrec.contents.nstxtcout = c_int(0)
    bAccept = bool(libmdrun.copy_and_add_molecule(byref(currmc),
                                                  byref(trialmc),
                                                  c_int(Ncurr),
                                                  c_int(Ntrial),
                                                  xadd,
                                                  vadd,
                                                  c_atomtype,
                                                  c_int(NAtomsOfMolType),
                                                  c_int(IndexOffSet),
                                                  SigmaCrossTerms2Frac,
                                                  c_int(len(SigmasAtomCrossTerms2Frac))))

    return bAccept

def MCInsertN0(trialmc = gt.mdrunner_comm,
               Ntrial = int,
               NAtomsOfMolType = int,
               SigmaVel = [],
               SigmasAtomCrossTerms2Frac = [[]],
               IndexOffSet = int,
               AtomType = [],
               Conf = []):

    c_intarray = c_int*NAtomsOfMolType
    c_atomtype = c_intarray()
    for i in range(NAtomsOfMolType):
        c_atomtype[i] = c_int(AtomType[i])

    c_realarray = c_real*len(SigmasAtomCrossTerms2Frac)*len(SigmasAtomCrossTerms2Frac)
    SigmaCrossTerms2Frac = c_realarray()
    for i in range(len(SigmasAtomCrossTerms2Frac)):
        for j in range(len(SigmasAtomCrossTerms2Frac)):
            SigmaCrossTerms2Frac[i][j] = c_real(SigmasAtomCrossTerms2Frac[i][j])

    rvecarray = rvec*NAtomsOfMolType
    vecx      = rvecarray()
    vecv      = rvecarray()
    r         = [0.0]*3
    r[gt.XX]  = random.random()*trialmc.box[gt.XX][gt.XX]
    r[gt.YY]  = random.random()*trialmc.box[gt.YY][gt.YY]
    r[gt.ZZ]  = random.random()*trialmc.box[gt.ZZ][gt.ZZ]
    ConfNew   = []
    if(NAtomsOfMolType>1):
        ConfNew = RotateMoleculeRandomly(Conf)
    else:
        ConfNew = Conf

#   maybe COM vel + random rotations?
    for i in range(NAtomsOfMolType):
        vecx[i][gt.XX] = r[gt.XX] + ConfNew[i][gt.XX]
        vecx[i][gt.YY] = r[gt.YY] + ConfNew[i][gt.YY]
        vecx[i][gt.ZZ] = r[gt.ZZ] + ConfNew[i][gt.ZZ]
        vecv[i][gt.XX] = random.gauss(0.0,SigmaVel[i])
        vecv[i][gt.YY] = random.gauss(0.0,SigmaVel[i])
        vecv[i][gt.ZZ] = random.gauss(0.0,SigmaVel[i])

    xadd = rvecarray()
    vadd = rvecarray()
    for i in range(NAtomsOfMolType):
        xadd[i][gt.XX] = c_real(vecx[i][gt.XX])
        xadd[i][gt.YY] = c_real(vecx[i][gt.YY])
        xadd[i][gt.ZZ] = c_real(vecx[i][gt.ZZ])
        vadd[i][gt.XX] = c_real(vecv[i][gt.XX])
        vadd[i][gt.YY] = c_real(vecv[i][gt.YY])
        vadd[i][gt.ZZ] = c_real(vecv[i][gt.ZZ])

    trialmc.inputrec.contents.nsteps    = c_int(0)
    trialmc.inputrec.contents.nstlog    = c_int(0)
    trialmc.inputrec.contents.nstxout   = c_int(0)
    trialmc.inputrec.contents.nstvout   = c_int(0)
    trialmc.inputrec.contents.nstfout   = c_int(0)
    trialmc.inputrec.contents.nstenergy = c_int(0)
    trialmc.inputrec.contents.nstxtcout = c_int(0)

    for i in range(NAtomsOfMolType):
        trialmc.state.contents.x[i][gt.XX] = xadd[i][gt.XX]
        trialmc.state.contents.x[i][gt.YY] = xadd[i][gt.YY]
        trialmc.state.contents.x[i][gt.ZZ] = xadd[i][gt.ZZ]
        trialmc.state.contents.v[i][gt.XX] = vadd[i][gt.XX]
        trialmc.state.contents.v[i][gt.YY] = vadd[i][gt.YY]
        trialmc.state.contents.v[i][gt.ZZ] = vadd[i][gt.ZZ]

    bAccept = True

    return bAccept

def MCRemove(currmc=gt.mdrunner_comm,
             trialmc=gt.mdrunner_comm,
             Ncurr=int,
             Ntrial=int,
             bGenerate=bool,
             iRemoveInput=int,
             IndexOffSet=int,
             NAtomsOfMolType=int):
    trialmc.inputrec.contents.nsteps    = c_int(0)
    trialmc.inputrec.contents.nstlog    = c_int(0)
    trialmc.inputrec.contents.nstxout   = c_int(0)
    trialmc.inputrec.contents.nstvout   = c_int(0)
    trialmc.inputrec.contents.nstfout   = c_int(0)
    trialmc.inputrec.contents.nstenergy = c_int(0)
    trialmc.inputrec.contents.nstxtcout = c_int(0)
#   select particle at random
    if(bGenerate):
        iRemove = random.randint(0,Ncurr-1)
    else:
        iRemove = iRemoveInput
    iRemove = 239
    libmdrun.copy_and_remove_molecule(byref(currmc),
                                      byref(trialmc),
                                      c_int(Ncurr),
                                      c_int(Ntrial),
                                      c_int(iRemove),
                                      c_int(IndexOffSet),
                                      c_int(NAtomsOfMolType))
    return iRemove

def GetVolRed(mc = gt.mdrunner_comm,Sigma = float):
    Vol  = mc.box[gt.XX][gt.XX]/Sigma
    Vol *= mc.box[gt.YY][gt.YY]/Sigma
    Vol *= mc.box[gt.ZZ][gt.ZZ]/Sigma
    return Vol

def GetVol(mc = gt.mdrunner_comm):
    Vol  = mc.box[gt.XX][gt.XX]
    Vol *= mc.box[gt.YY][gt.YY]
    Vol *= mc.box[gt.ZZ][gt.ZZ]
    return Vol

def GetIndexForPotentialEnergy(mc = gt.mdrunner_comm):
    nener = mc.mdebin.contents.ebin.contents.nener
    for i in range(nener):
        if(mc.mdebin.contents.ebin.contents.enm[i]=='Potential'):
            iPot = i
    return iPot

def GetUPotFromMdebin(mc = gt.mdrunner_comm,iPot = int):
    return float(mc.mdebin.contents.ebin.contents.e[iPot].e)

def GetUPotOfTrialRemoval(mc = gt.mdrunner_comm):
    nener = mc.mdebin.contents.ebin.contents.nener
    UPot  = 0.0
    for i in range(nener):
        if(re.search('WTrial-rest',mc.mdebin.contents.ebin.contents.enm[i])):
            UPot += float(mc.mdebin.contents.ebin.contents.e[i].e)
    return UPot

def GetLJC6(mc = gt.mdrunner_comm):
    return mc.mtop.contents.ffparams.iparams.contents.lj.c6

def GetEps():
    return 5.0 # MARTINI eps of water (kJ/mol)

def GetSigma(mc = gt.mdrunner_comm):
    return round(math.pow(GetLJC6(mc)/(4.0*GetEps()),(1.0/6.0)),2)

def GetSigmas(Mc = gt.mdrunner_comm):
    AtomSigmas = []
    for i in range(Mc.mtop.contents.ffparams.atnr):
        LjC6  = Mc.mtop.contents.ffparams.iparams[i].lj.c6
        LjC12 = Mc.mtop.contents.ffparams.iparams[i].lj.c12
        if(LjC6==0.0 and LjC12==0.0):
            AtomSigmas.append(0.0)
        else:
            AtomSigmas.append(math.pow(LjC12/LjC6,(1.0/6.0)))

    return AtomSigmas

def GetSigmaCrossTerms(AtomSigmas=[],Fraction=float):
#    AtomSigmaCrossTerms  = [[]*NEntries]*NEntries
#    AtomSigmaCrossTerms2 = [[]*NEntries]*NEntries
    AtomSigmaCrossTerms      = []
    AtomSigmaCrossTerms2     = []
    AtomSigmaCrossTerms2Frac = []
    for i in range(len(AtomSigmas)):
        List      = []
        List2     = []
        List2Frac = []
        for j in range(len(AtomSigmas)):
            CrossTerm = (AtomSigmas[i]+AtomSigmas[j])*0.5
            List.append(CrossTerm)
            List2.append(math.pow(CrossTerm,2.0))
            List2Frac.append(math.pow(CrossTerm*Fraction,2.0))
        AtomSigmaCrossTerms.append(List)
        AtomSigmaCrossTerms2.append(List2)
        AtomSigmaCrossTerms2Frac.append(List2Frac)
        del List,List2,List2Frac

    return AtomSigmaCrossTerms,AtomSigmaCrossTerms2,AtomSigmaCrossTerms2Frac


#4 eps ((sig)^12/(r)^12-(sig)^6/(r)^6)
#c12 = 4 eps sig^12
#c6  = 4 eps sig^6

#    for i in range(Mc.mtop.contents.nmoltype):
#        SigmaList = []
#        for j in in range(Mc.mtop.contents.moltype[i].atoms.nr)
#            print i,j,SigmaList.append()
#        AtomSigmas.append(SigmaList)
#        del SigmaList

    return AtomSigmas

def RemoveBaks(pwd):
    for file in os.listdir(pwd):
        if(file[0]=='#' and file[-1]=='#'):
            os.remove(pwd+'/'+file)

def RenameOutputFiles(pwd):
    FileList = ['ener.edr','traj.trr','traj.xtc','counfout.gro']
    for File in FileList:
        CompletePath = pwd+'/'+File
        if(os.path.isfile(CompletePath)):
            os.rename(CompletePath, CompletePath+'.bak')

def InitMdrunnerComm():
    mc             = gt.mdrunner_comm()
    mc.state       = gt.NULL
    mc.buf         = gt.NULL
    mc.f           = gt.NULL
    mc.mtop        = gt.NULL
    mc.mdatoms     = gt.NULL
    mc.fr          = gt.NULL
    mc.fcd         = gt.NULL
    mc.ewaldcoeff  = c_real(0)
    mc.pmedata     = gt.NULL
    mc.vsite       = gt.NULL
    mc.nsteps_done = c_int(0)
    mc.mdebin      = gt.NULL
    mc.nsteps_done = c_int(0)
    mc.xsave       = gt.NULL
    mc.vsave       = gt.NULL
    return mc

def GetNMolOfMolType(Mc = gt.mdrunner_comm, GCMCMolType = str):
    GCMCMolTypeIndex = 0
    for i in range(Mc.mtop.contents.nmoltype):
        if(Mc.mtop.contents.moltype[i].name.contents.value==GCMCMolType):
            GCMCMolTypeIndex = i
            break
    GCMCMolBlockIndex = 0
    for i in range(Mc.mtop.contents.nmolblock):
        if(Mc.mtop.contents.molblock[i].type==GCMCMolTypeIndex):
            GCMCMolBlockIndex = i
            break
    return Mc.mtop.contents.molblock[GCMCMolBlockIndex].nmol

def GetMassOfMolType(Mc = gt.mdrunner_comm, GCMCMolType = str):
    GCMCMolTypeIndex = 0
    for i in range(Mc.mtop.contents.nmoltype):
        if(Mc.mtop.contents.moltype[i].name.contents.value==GCMCMolType):
            GCMCMolTypeIndex = i
            break
    NAtomsOfMolType = Mc.mtop.contents.moltype[GCMCMolTypeIndex].atoms.nr
    Mass     = 0.0
    AtomMass = []
    for i in range(NAtomsOfMolType):
        Mass += Mc.mtop.contents.moltype[GCMCMolTypeIndex].atoms.atom[i].m
        AtomMass.append(Mc.mtop.contents.moltype[GCMCMolTypeIndex].atoms.atom[i].m)

    return Mass,AtomMass,NAtomsOfMolType

def MCInit(mc = gt.mdrunner_comm, GCMCMolType = str):
    GCMCMolTypeIndex = 0
    for i in range(mc.mtop.contents.nmoltype):
        if(mc.mtop.contents.moltype[i].name.contents.value==GCMCMolType):
            GCMCMolTypeIndex = i
            break
    GCMCMolBlockIndex = 0
    for i in range(mc.mtop.contents.nmolblock):
        if(mc.mtop.contents.molblock[i].type==GCMCMolTypeIndex):
            GCMCMolBlockIndex = i
            break
    return GCMCMolTypeIndex,GCMCMolBlockIndex

# from mdrun.h
MD_GLAS         = (1<<1)
MD_IONIZE       = (1<<3)
MD_RERUN        = (1<<4)
MD_RERUN_VSITE  = (1<<5)
MD_SEPPOT       = (1<<7)
MD_PARTDEC      = (1<<9)
MD_DDBONDCHECK  = (1<<10)
MD_DDBONDCOMM   = (1<<11)
MD_CONFOUT      = (1<<12)
MD_NOGSTAT      = (1<<13)
MD_REPRODUCIBLE = (1<<14)
MD_READ_RNG     = (1<<15)
MD_APPENDFILES  = (1<<16)
MD_READ_EKIN    = (1<<17)
MD_STARTFROMCPT = (1<<18)

# from mdrun.h (octals)
LIST_SCALARS  = 0001
LIST_INPUTREC = 0002
LIST_TOP      = 0004
LIST_X        = 0010
LIST_V        = 0020
LIST_F        = 0040
LIST_LOAD     = 0100

def hybrid(Mu,NStart,NMax,NMin,TprDir,MolName):
    MuRed = Mu*iEps
    (ddnoSEL,ddnoINTERLEAVE,ddnoPP_PME,ddnoCARTESIAN,ddnoNR) = map(c_int,xrange(5))

    descarray = c_char_p*191
    desc      = descarray(
        c_char_p("The mdrun program is the main computational chemistry engine"),
        c_char_p("within GROMACS. Obviously, it performs Molecular Dynamics simulations,"),
        c_char_p("but it can also perform Stochastic Dynamics, Energy Minimization,"),
        c_char_p("test particle insertion or (re)calculation of energies."),
        c_char_p("Normal mode analysis is another option. In this case mdrun"),
        c_char_p("builds a Hessian matrix from single conformation."),
        c_char_p("For usual Normal Modes-like calculations, make sure that"),
        c_char_p("the structure provided is properly energy-minimized."),
        c_char_p("The generated matrix can be diagonalized by g_nmeig.[PAR]"),
        c_char_p("The mdrun program reads the run input file ([TT]-s[tt])"),
        c_char_p("and distributes the topology over nodes if needed."),
        c_char_p("mdrun produces at least four output files."),
        c_char_p("A single log file ([TT]-g[tt]) is written, unless the option"),
        c_char_p("[TT]-seppot[tt] is used, in which case each node writes a log file."),
        c_char_p("The trajectory file ([TT]-o[tt]), contains coordinates, velocities and"),
        c_char_p("optionally forces."),
        c_char_p("The structure file ([TT]-c[tt]) contains the coordinates and"),
        c_char_p("velocities of the last step."),
        c_char_p("The energy file ([TT]-e[tt]) contains energies, the temperature,"),
        c_char_p("pressure, etc, a lot of these things are also printed in the log file."),
        c_char_p("Optionally coordinates can be written to a compressed trajectory file"),
        c_char_p("([TT]-x[tt]).[PAR]"),
        c_char_p("The option [TT]-dgdl[tt] is only used when free energy perturbation is"),
        c_char_p("turned on.[PAR]"),
        c_char_p("When mdrun is started using MPI with more than 1 node, parallelization"),
        c_char_p("is used. By default domain decomposition is used, unless the [TT]-pd[tt]"),
        c_char_p("option is set, which selects particle decomposition.[PAR]"),
        c_char_p("With domain decomposition, the spatial decomposition can be set"),
        c_char_p("with option [TT]-dd[tt]. By default mdrun selects a good decomposition."),
        c_char_p("The user only needs to change this when the system is very inhomogeneous."),
        c_char_p("Dynamic load balancing is set with the option [TT]-dlb[tt],"),
        c_char_p("which can give a significant performance improvement,"),
        c_char_p("especially for inhomogeneous systems. The only disadvantage of"),
        c_char_p("dynamic load balancing is that runs are no longer binary reproducible,"),
        c_char_p("but in most cases this is not important."),
        c_char_p("By default the dynamic load balancing is automatically turned on"),
        c_char_p("when the measured performance loss due to load imbalance is 5% or more."),
        c_char_p("At low parallelization these are the only important options"),
        c_char_p("for domain decomposition."),
        c_char_p("At high parallelization the options in the next two sections"),
        c_char_p("could be important for increasing the performace."),
        c_char_p("[PAR]"),
        c_char_p("When PME is used with domain decomposition, separate nodes can"),
        c_char_p("be assigned to do only the PME mesh calculation;"),
        c_char_p("this is computationally more efficient starting at about 12 nodes."),
        c_char_p("The number of PME nodes is set with option [TT]-npme[tt],"),
        c_char_p("this can not be more than half of the nodes."),
        c_char_p("By default mdrun makes a guess for the number of PME"),
        c_char_p("nodes when the number of nodes is larger than 11 or performance wise"),
        c_char_p("not compatible with the PME grid x dimension."),
        c_char_p("But the user should optimize npme. Performance statistics on this issue"),
        c_char_p("are written at the end of the log file."),
        c_char_p("For good load balancing at high parallelization, the PME grid x and y"),
        c_char_p("dimensions should be divisible by the number of PME nodes"),
        c_char_p("(the simulation will run correctly also when this is not the case)."),
        c_char_p("[PAR]"),
        c_char_p("This section lists all options that affect the domain decomposition."),
        c_char_p("[BR]"),
        c_char_p("Option [TT]-rdd[tt] can be used to set the required maximum distance"),
        c_char_p("for inter charge-group bonded interactions."),
        c_char_p("Communication for two-body bonded interactions below the non-bonded"),
        c_char_p("cut-off distance always comes for free with the non-bonded communication."),
        c_char_p("Atoms beyond the non-bonded cut-off are only communicated when they have"),
        c_char_p("missing bonded interactions; this means that the extra cost is minor"),
        c_char_p("and nearly indepedent of the value of [TT]-rdd[tt]."),
        c_char_p("With dynamic load balancing option [TT]-rdd[tt] also sets"),
        c_char_p("the lower limit for the domain decomposition cell sizes."),
        c_char_p("By default [TT]-rdd[tt] is determined by mdrun based on"),
        c_char_p("the initial coordinates. The chosen value will be a balance"),
        c_char_p("between interaction range and communication cost."),
        c_char_p("[BR]"),
        c_char_p("When inter charge-group bonded interactions are beyond"),
        c_char_p("the bonded cut-off distance, mdrun terminates with an error message."),
        c_char_p("For pair interactions and tabulated bonds"),
        c_char_p("that do not generate exclusions, this check can be turned off"),
        c_char_p("with the option [TT]-noddcheck[tt]."),
        c_char_p("[BR]"),
        c_char_p("When constraints are present, option [TT]-rcon[tt] influences"),
        c_char_p("the cell size limit as well."),
        c_char_p("Atoms connected by NC constraints, where NC is the LINCS order plus 1,"),
        c_char_p("should not be beyond the smallest cell size. A error message is"),
        c_char_p("generated when this happens and the user should change the decomposition"),
        c_char_p("or decrease the LINCS order and increase the number of LINCS iterations."),
        c_char_p("By default mdrun estimates the minimum cell size required for P-LINCS"),
        c_char_p("in a conservative fashion. For high parallelization it can be useful"),
        c_char_p("to set the distance required for P-LINCS with the option [TT]-rcon[tt]."),
        c_char_p("[BR]"),
        c_char_p("The [TT]-dds[tt] option sets the minimum allowed x, y and/or z scaling"),
        c_char_p("of the cells with dynamic load balancing. mdrun will ensure that"),
        c_char_p("the cells can scale down by at least this factor. This option is used"),
        c_char_p("for the automated spatial decomposition (when not using [TT]-dd[tt])"),
        c_char_p("as well as for determining the number of grid pulses, which in turn"),
        c_char_p("sets the minimum allowed cell size. Under certain circumstances"),
        c_char_p("the value of [TT]-dds[tt] might need to be adjusted to account for"),
        c_char_p("high or low spatial inhomogeneity of the system."),
        c_char_p("[PAR]"),
        c_char_p("The option [TT]-nosum[tt] can be used to only sum the energies"),
        c_char_p("at every neighbor search step and energy output step."),
        c_char_p("This can improve performance for highly parallel simulations"),
        c_char_p("where this global communication step becomes the bottleneck."),
        c_char_p("For a global thermostat and/or barostat the temperature"),
        c_char_p("and/or pressure will also only be updated every nstlist steps."),
        c_char_p("With this option the energy file will not contain averages and"),
        c_char_p("fluctuations over all integration steps.[PAR]"),
        c_char_p("With [TT]-rerun[tt] an input trajectory can be given for which "),
        c_char_p("forces and energies will be (re)calculated. Neighbor searching will be"),
        c_char_p("performed for every frame, unless [TT]nstlist[tt] is zero"),
        c_char_p("(see the [TT].mdp[tt] file).[PAR]"),
        c_char_p("ED (essential dynamics) sampling is switched on by using the [TT]-ei[tt]"),
        c_char_p("flag followed by an [TT].edi[tt] file."),
        c_char_p("The [TT].edi[tt] file can be produced using options in the essdyn"),
        c_char_p("menu of the WHAT IF program. mdrun produces a [TT].edo[tt] file that"),
        c_char_p("contains projections of positions, velocities and forces onto selected"),
        c_char_p("eigenvectors.[PAR]"),
        c_char_p("When user-defined potential functions have been selected in the"),
        c_char_p("[TT].mdp[tt] file the [TT]-table[tt] option is used to pass mdrun"),
        c_char_p("a formatted table with potential functions. The file is read from"),
        c_char_p("either the current directory or from the GMXLIB directory."),
        c_char_p("A number of pre-formatted tables are presented in the GMXLIB dir,"),
        c_char_p("for 6-8, 6-9, 6-10, 6-11, 6-12 Lennard Jones potentials with"),
        c_char_p("normal Coulomb."),
        c_char_p("When pair interactions are present a separate table for pair interaction"),
        c_char_p("functions is read using the [TT]-tablep[tt] option.[PAR]"),
        c_char_p("When tabulated bonded functions are present in the topology,"),
        c_char_p("interaction functions are read using the [TT]-tableb[tt] option."),
        c_char_p("For each different tabulated interaction type the table file name is"),
        c_char_p("modified in a different way: before the file extension an underscore is"),
        c_char_p("appended, then a b for bonds, an a for angles or a d for dihedrals"),
        c_char_p("and finally the table number of the interaction type.[PAR]"),
        c_char_p("The options [TT]-px[tt] and [TT]-pf[tt] are used for writing pull COM"),
        c_char_p("coordinates and forces when pulling is selected"),
        c_char_p("in the [TT].mdp[tt] file.[PAR]"),
        c_char_p("With [TT]-multi[tt] multiple systems are simulated in parallel."),
        c_char_p("As many input files are required as the number of systems."),
        c_char_p("The system number is appended to the run input and each output filename,"),
        c_char_p("for instance topol.tpr becomes topol0.tpr, topol1.tpr etc."),
        c_char_p("The number of nodes per system is the total number of nodes"),
        c_char_p("divided by the number of systems."),
        c_char_p("One use of this option is for NMR refinement: when distance"),
        c_char_p("or orientation restraints are present these can be ensemble averaged"),
        c_char_p("over all the systems.[PAR]"),
        c_char_p("With [TT]-replex[tt] replica exchange is attempted every given number"),
        c_char_p("of steps. The number of replicas is set with the [TT]-multi[tt] option,"),
        c_char_p("see above."),
        c_char_p("All run input files should use a different coupling temperature,"),
        c_char_p("the order of the files is not important. The random seed is set with"),
        c_char_p("[TT]-reseed[tt]. The velocities are scaled and neighbor searching"),
        c_char_p("is performed after every exchange.[PAR]"),
        c_char_p("Finally some experimental algorithms can be tested when the"),
        c_char_p("appropriate options have been given. Currently under"),
        c_char_p("investigation are: polarizability, glass simulations"),
        c_char_p("and X-Ray bombardments."),
        c_char_p("[PAR]"),
        c_char_p("The option [TT]-pforce[tt] is useful when you suspect a simulation"),
        c_char_p("crashes due to too large forces. With this option coordinates and"),
        c_char_p("forces of atoms with a force larger than a certain value will"),
        c_char_p("be printed to stderr."),
        c_char_p("[PAR]"),
        c_char_p("Checkpoints containing the complete state of the system are written"),
        c_char_p("at regular intervals (option [TT]-cpt[tt]) to the file [TT]-cpo[tt],"),
        c_char_p("unless option [TT]-cpt[tt] is set to -1."),
        c_char_p("A simulation can be continued by reading the full state from file"),
        c_char_p("with option [TT]-cpi[tt]. This option is intelligent in the way that"),
        c_char_p("if no checkpoint file is found, Gromacs just assumes a normal run and"),
        c_char_p("starts from the first step of the tpr file."),
        c_char_p("The simulation part number is added to all output files,"),
        c_char_p("unless [TT]-append[tt] or [TT]-noaddpart[tt] are set."),
        c_char_p("[PAR]"),
        c_char_p("With checkpointing you can also use the option [TT]-append[tt] to"),
        c_char_p("just continue writing to the previous output files. This is not"),
        c_char_p("enabled by default since it is potentially dangerous if you move files,"),
        c_char_p("but if you just leave all your files in place and restart mdrun with"),
        c_char_p("exactly the same command (with options [TT]-cpi[tt] and [TT]-append[tt])"),
        c_char_p("the result will be the same as from a single run. The contents will"),
        c_char_p("be binary identical (unless you use dynamic load balancing),"),
        c_char_p("but for technical reasons there might be some extra energy frames when"),
        c_char_p("using checkpointing (necessary for restarts without appending)."),
        c_char_p("[PAR]"),
        c_char_p("With option [TT]-maxh[tt] a simulation is terminated and a checkpoint"),
        c_char_p("file is written at the first neighbor search step where the run time"),
        c_char_p("exceeds [TT]-maxh[tt]*0.99 hours."),
        c_char_p("[PAR]"),
        c_char_p("When mdrun receives a TERM signal, it will set nsteps to the current"),
        c_char_p("step plus one. When mdrun receives a USR1 signal, it will stop after"),
        c_char_p("the next neighbor search step (with nstlist=0 at the next step)."),
        c_char_p("In both cases all the usual output will be written to file."),
        c_char_p("When running with MPI, a signal to one of the mdrun processes"),
        c_char_p("is sufficient, this signal should not be sent to mpirun or"),
        c_char_p("the mdrun process that is the parent of the others."),
        c_char_p("[PAR]"),
        c_char_p("When mdrun is started with MPI, it does not run niced by default.")
        )

    pt_commrec = POINTER(gt.t_commrec)
    cr         = pt_commrec

    fnmarray = gt.t_filenm * 27
    fnm      = fnmarray(gt.t_filenm(ftp=gt.efTPX, opt=c_char_p(),          fn=c_char_p(),           flag=gt.ffREAD , nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efTRN, opt=c_char_p("-o"),      fn=c_char_p(),           flag=gt.ffWRITE, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXTC, opt=c_char_p("-x"),      fn=c_char_p(),           flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efCPT, opt=c_char_p("-cpi"),    fn=c_char_p(),           flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efCPT, opt=c_char_p("-cpo"),    fn=c_char_p(),           flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efSTO, opt=c_char_p("-c"),      fn=c_char_p("confout"),  flag=gt.ffWRITE, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efENX, opt=c_char_p("-e"),      fn=c_char_p("ener"),     flag=gt.ffWRITE, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efLOG, opt=c_char_p("-g"),      fn=c_char_p("md"),       flag=gt.ffWRITE, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-dgdl"),   fn=c_char_p("dgdl"),     flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-field"),  fn=c_char_p("field"),    flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-table"),  fn=c_char_p("table"),    flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-tablep"), fn=c_char_p("tablep"),   flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-tableb"), fn=c_char_p("table"),    flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efTRX, opt=c_char_p("-rerun"),  fn=c_char_p("rerun"),    flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-tpi"),    fn=c_char_p("tpi"),      flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-tpid"),   fn=c_char_p("tpidist"),  flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efEDI, opt=c_char_p("-ei"),     fn=c_char_p("sam"),      flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efEDO, opt=c_char_p("-eo"),     fn=c_char_p("sam"),      flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efGCT, opt=c_char_p("-j"),      fn=c_char_p("wham"),     flag=gt.ffOPTRD, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efGCT, opt=c_char_p("-jo"),     fn=c_char_p("bam"),      flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-ffout"),  fn=c_char_p("gct"),      flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-devout"), fn=c_char_p("deviatie"), flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-runav"),  fn=c_char_p("runaver"),  flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-px"),     fn=c_char_p("pullx"),    flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efXVG, opt=c_char_p("-pf"),     fn=c_char_p("pullf"),    flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efMTX, opt=c_char_p("-mtx"),    fn=c_char_p("nm"),       flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()),
                        gt.t_filenm(ftp=gt.efNDX, opt=c_char_p("-dn"),     fn=c_char_p("dipole"),   flag=gt.ffOPTWR, nfiles=c_int(0), fns=POINTER(c_char_p)()))
    NFILE =  c_int(len(fnm))

    bCart         = c_int(gt.FALSE)
    bPPPME        = c_int(gt.FALSE)
    bPartDec      = c_int(gt.FALSE)
    bDDBondCheck  = c_int(gt.TRUE)
    bDDBondComm   = c_int(gt.TRUE)
    bSumEner      = c_int(gt.TRUE)
    bVerbose      = c_int(gt.FALSE)
    bCompact      = c_int(gt.TRUE)
    bSepPot       = c_int(gt.FALSE)
    bRerunVSite   = c_int(gt.FALSE)
    bGlas         = c_int(gt.FALSE)
    bIonize       = c_int(gt.FALSE)
    bConfout      = c_int(gt.TRUE)
    bReproducible = c_int(gt.FALSE)

    npme         = c_int(-1)
    nmultisim    = c_int(0)
    repl_ex_nst  = c_int(0)
    repl_ex_seed = c_int(-1)
    nstepout     = c_int(100)
    nthreads     = c_int(1)

    realddxyz  = rvec(c_real(0.0),c_real(0.0),c_real(0.0))
    ddnoarray  = c_char_p*(ddnoNR.value+1)
    ddno_opt   = ddnoarray(c_char_p(gt.NULL),
                           c_char_p("interleave"),
                           c_char_p("pp_pme"),
                           c_char_p("cartesian"),
                           c_char_p(gt.NULL))
    ddndlbarray = c_char_p*5
    dddlb_opt   = ddndlbarray(c_char_p(gt.NULL),
                              c_char_p("auto"),
                              c_char_p("no"),
                              c_char_p("yes"),
                              c_char_p(gt.NULL))

    rdd          = c_real(0.0)
    rconstr      = c_real(0.0)
    dlb_scale    = c_real(0.8)
    pforce       = c_real(-1)
    ddcsx        = c_char_p(gt.NULL)
    ddcsy        = c_char_p(gt.NULL)
    ddcsz        = c_char_p(gt.NULL)
    cpt_period   = c_real(15.0)
    max_hours    = c_real(-1)
    bAppendFiles = c_int(gt.FALSE)
    bAddPart     = c_int(gt.TRUE)
    t_pargsarray = gt.t_pargs*32
    pa           = t_pargsarray(
                        gt.t_pargs(option=c_char_p("-pd"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bPartDec)),
                                   desc=c_char_p("Use particle decompostion")),
                        gt.t_pargs(option=c_char_p("-dd"),
                                   bSet=gt.FALSE,
                                   type=gt.etRVEC,
                                   u=gt.u_u(addressof(realddxyz)),
                                   desc=c_char_p("Domain decomposition grid, 0 is optimize")),
                        gt.t_pargs(option=c_char_p("-nt"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(nthreads)),
                                   desc=c_char_p("HIDDENNumber of threads to start on each node")),
                        gt.t_pargs(option=c_char_p("-npme"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(npme)),
                                   desc=c_char_p("Number of separate nodes to be used for PME, -1 is guess")),
                        gt.t_pargs(option=c_char_p("-ddorder"),
                                   bSet=gt.FALSE,
                                   type=gt.etENUM,
                                   u=gt.u_u(addressof(ddno_opt)),
                                   desc=c_char_p("DD node order")),
                        gt.t_pargs(option=c_char_p("-ddcheck"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bDDBondCheck)),
                                   desc=c_char_p("Check for all bonded interactions with DD")),
                        gt.t_pargs(option=c_char_p("-ddbondcomm"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bDDBondComm)),
                                   desc=c_char_p("HIDDENUse special bonded atom communication when -rdd > cut-off")),
                        gt.t_pargs(option=c_char_p("-rdd"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(rdd)),
                                   desc=c_char_p("The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates")),
                        gt.t_pargs(option=c_char_p("-rcon"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(rconstr)),
                                   desc=c_char_p("Maximum distance for P-LINCS (nm), 0 is estimate")),
                        gt.t_pargs(option=c_char_p("-dlb"),
                                   bSet=gt.FALSE,
                                   type=gt.etENUM,
                                   u=gt.u_u(addressof(dddlb_opt)),
                                   desc=c_char_p("Dynamic load balancing (with DD)")),
                        gt.t_pargs(option=c_char_p("-dds"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(dlb_scale)),
                                   desc=c_char_p("Minimum allowed dlb scaling of the DD cell size")),
                        gt.t_pargs(option=c_char_p("-ddcsx"),
                                   bSet=gt.FALSE,
                                   type=gt.etSTR,
                                   u=gt.u_u(addressof(ddcsx)),
                                   desc=c_char_p("HIDDENThe DD cell sizes in x")),
                        gt.t_pargs(option=c_char_p("-ddcsy"),
                                   bSet=gt.FALSE,
                                   type=gt.etSTR,
                                   u=gt.u_u(addressof(ddcsy)),
                                   desc=c_char_p("HIDDENThe DD cell sizes in y")),
                        gt.t_pargs(option=c_char_p("-ddcsz"),
                                   bSet=gt.FALSE,
                                   type=gt.etSTR,
                                   u=gt.u_u(addressof(ddcsz)),
                                   desc=c_char_p("HIDDENThe DD cell sizes in z")),
                        gt.t_pargs(option=c_char_p("-sum"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bSumEner)),
                                   desc=c_char_p("Sum the energies at every step")),
                        gt.t_pargs(option=c_char_p("-v"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bVerbose)),
                                   desc=c_char_p("Be loud and noisy")),
                        gt.t_pargs(option=c_char_p("-compact"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bCompact)),
                                   desc=c_char_p("Write a compact log file")),
                        gt.t_pargs(option=c_char_p("-seppot"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bSepPot)),
                                   desc=c_char_p("Write separate V and dVdl terms for each interaction type and node to the log file(s)")),
                        gt.t_pargs(option=c_char_p("-pforce"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(pforce)),
                                   desc=c_char_p("Print all forces larger than this (kJ/mol nm)")),
                        gt.t_pargs(option=c_char_p("-reprod"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bReproducible)),
                                   desc=c_char_p("Try to avoid optimizations that affect binary reproducibility")),
                        gt.t_pargs(option=c_char_p("-cpt"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(cpt_period)),
                                   desc=c_char_p("Checkpoint interval (minutes)")),
                        gt.t_pargs(option=c_char_p("-append"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bAppendFiles)),
                                   desc=c_char_p("Append to previous output files when continuing from checkpoint")),
                        gt.t_pargs(option=c_char_p("-addpart"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bAddPart)),
                                   desc=c_char_p("Add the simulation part number to all output files when continuing from checkpoint")),
                        gt.t_pargs(option=c_char_p("-maxh"),
                                   bSet=gt.FALSE,
                                   type=gt.etREAL,
                                   u=gt.u_u(addressof(max_hours)),
                                   desc=c_char_p("Terminate after 0.99 times this time (hours)")),
                        gt.t_pargs(option=c_char_p("-multi"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(nmultisim)),
                                   desc=c_char_p("Do multiple simulations in parallel")),
                        gt.t_pargs(option=c_char_p("-replex"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(repl_ex_nst)),
                                   desc=c_char_p("Attempt replica exchange every # steps")),
                        gt.t_pargs(option=c_char_p("-reseed"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(repl_ex_seed)),
                                   desc=c_char_p("Seed for replica exchange, -1 is generate a seed")),
                        gt.t_pargs(option=c_char_p("-rerunvsite"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bRerunVSite)),
                                   desc=c_char_p("HIDDENRecalculate virtual site coordinates with -rerun")),
                        gt.t_pargs(option=c_char_p("-glas"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bGlas)),
                                   desc=c_char_p("Do glass simulation with special long range corrections")),
                        gt.t_pargs(option=c_char_p("-ionize"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bIonize)),
                                   desc=c_char_p("Do a simulation including the effect of an X-Ray bombardment on your system")),
                        gt.t_pargs(option=c_char_p("-confout"),
                                   bSet=gt.FALSE,
                                   type=gt.etBOOL,
                                   u=gt.u_u(addressof(bConfout)),
                                   desc=c_char_p("HIDDENWrite the last configuration with -c")),
                        gt.t_pargs(option=c_char_p("-stepout"),
                                   bSet=gt.FALSE,
                                   type=gt.etINT,
                                   u=gt.u_u(addressof(nstepout)),
                                   desc=c_char_p("HIDDENFrequency of writing the remaining runtime"))
                        )

    ed             = gt.gmx_edsam_t()
    Flags          = c_ulong()
    PCA_Flags      = c_ulong()
    ddxyz          = ivec()
    dd_node_order  = c_int()
    HaveCheckpoint = c_int()
    fplog          = gt.FILE_p()
    fptes          = gt.FILE_p()
    sim_part       = c_int()
    suffixarray    = c_char*gt.STRLEN
    suffix         = suffixarray()

############################
#   Quick hack: No options #
############################
    argv_array = c_char_p*1
    argv       = argv_array(c_char_p("mdrun"))
    argc       = c_int(1)

    init_par         = libmdrun.init_par
    init_par.restype = pt_commrec
    cr               = (init_par(byref(argc),byref(argv)))

    if(commrec.MASTER(cr)):
        libmdrun.CopyRight(stderr,argv[0])

    PCA_Flags = (statutil.PCA_KEEP_ARGS | statutil.PCA_NOEXIT_ON_ARGS | statutil.PCA_CAN_SET_DEFFNM | (0 if commrec.MASTER(cr) else statutil.PCA_QUIET))
    if (not c_int.in_dll(libmdrun,"gmx_parallel_env").value):
        PCA_Flags = (PCA_Flags | statutil.PCA_BE_NICE)
    PCA_Flags = c_ulong(PCA_Flags)

    libmdrun.parse_common_args(byref(argc),\
                               argv,\
                               PCA_Flags,\
                               NFILE,\
                               fnm,\
                               c_int(len(pa)),\
                               pa,\
                               c_int(len(desc)),\
                               desc,\
                               c_int(0),\
                               POINTER(c_char_p)())

    dd_node_order = c_int(libmdrun.nenum(ddno_opt))
    cr.contents.npmenodes = npme

    if(nthreads.value>1):
        print "GROMACS compiled without threads support - can only use one thread"
        sys.exit(1)

    if(repl_ex_nst.value!=0 and nmultisim.value<2):
        print "Need at least two replicas for replica exchange (option -multi)"
        sys.exit(1)

    if(nmultisim.value>1 and commrec.PAR(cr)):
        libmdrun.init_multisystem(cr,nmultisim,NFILE,fnm,c_int(gt.TRUE))

#  /* Check if there is ANY checkpoint file available */
    sim_part = c_int(1)
    if(libmdrun.opt2bSet(c_char_p("-cpi"),NFILE,fnm)):
        sim_part = libmdrun.read_checkpoint_simulation_part(libmdrun.opt2fn("-cpi",NFILE,fnm)) + 1
#       /* sim_part will now be 1 if no checkpoint file was found */
        if(sim_part==1 and commrec.MASTER(cr)):
            print "No previous checkpoint file present, assuming this is a new run.\n"

    if(sim_part.value<=1):
        bAppendFiles = c_int(gt.FALSE)

    if((not bAppendFiles.value) and bAddPart.value and sim_part.value>1):
#        /* This is a continuation run, rename trajectory output files (except checkpoint files) */
#        /* create new part name first (zero-filled) */
        if(sim_part.value<10):
            libc.sprintf(suffix,"part000%d",sim_part)
        elif(sim_part.value<100):
            libc.sprintf(suffix,"part00%d",sim_part)
        elif(sim_part.value<1000):
            libc.sprintf(suffix,"part0%d",sim_part)
        else:
            libc.sprintf(suffix,"part%d",sim_part)
        libmdrun.add_suffix_to_output_names(fnm,NFILE,suffix)
        libc.fprintf(stdout,"Checkpoint file is from part %d, new output files will be suffixed %s.\n",c_int(sim_part.value-1),suffix)

    Flags = MD_RERUN if(libmdrun.opt2bSet(c_char_p("-rerun"),NFILE,fnm)) else 0
    Flags = Flags | (MD_SEPPOT if(bSepPot.value) else 0)
    Flags = Flags | (MD_IONIZE if(bIonize.value) else 0)
    Flags = Flags | (MD_GLAS if(bGlas.value) else 0)
    Flags = Flags | (MD_PARTDEC if(bPartDec.value) else 0)
    Flags = Flags | (MD_DDBONDCHECK if(bDDBondCheck.value) else 0)
    Flags = Flags | (MD_DDBONDCOMM if(bDDBondComm.value) else 0)
    Flags = Flags | (MD_CONFOUT if(bConfout.value) else 0)
    Flags = Flags | (MD_NOGSTAT if(not bSumEner.value) else 0)
    Flags = Flags | (MD_RERUN_VSITE if(bRerunVSite.value) else 0)
    Flags = Flags | (MD_REPRODUCIBLE if(bReproducible.value) else 0)
    Flags = Flags | (MD_APPENDFILES if(bAppendFiles.value) else 0)
    Flags = Flags | (MD_STARTFROMCPT if(sim_part.value>1) else 0)
    Flags = c_ulong(Flags)

#  /* We postpone opening the log file if we are appending, so we can first truncate
#   * the old log file and append to the correct position there instead.
#   */
    if(commrec.MASTER(cr) and (not bAppendFiles.value)):
        fplog = gt.FILE_p(libmdrun.gmx_log_open(c_char_p(libmdrun.ftp2fn(gt.efLOG,NFILE,fnm)),cr,c_int(not bSepPot.value),Flags))
        libmdrun.CopyRight(fplog,argv[0])
        libmdrun.please_cite(fplog,c_char_p("Hess2008b"))
        libmdrun.please_cite(fplog,c_char_p("Spoel2005a"))
        libmdrun.please_cite(fplog,c_char_p("Lindahl2001a"))
        libmdrun.please_cite(fplog,c_char_p("Berendsen95a"))
    else:
        fplog = gt.FILE_p()

#  /* Essential dynamics */
    if (libmdrun.opt2bSet(c_char_p("-ei"),NFILE,fnm)):
#   /* Open input and output files, allocate space for ED data structure */
        ed = libmdrun.ed_open(NFILE,fnm,cr)
    else:
        ed = gt.gmx_edsam_t()

    ddxyz[gt.XX] = c_int(int(realddxyz[gt.XX] + 0.5))
    ddxyz[gt.YY] = c_int(int(realddxyz[gt.YY] + 0.5))
    ddxyz[gt.ZZ] = c_int(int(realddxyz[gt.ZZ] + 0.5))

#####################
#    Initialization #
#####################
    os.system('ln -sf '+TprDir+'/W'+str(NStart)+'.tpr ./topol.tpr')
    CurrMc     = InitMdrunnerComm()
    DictMdrCmm = {}
    GCMCMolType = 'W' # input, for now single-atom particles only
    libmdrun.mdrunner_initialize(fplog,cr,NFILE,fnm,bVerbose,
                                 ddxyz,dd_node_order,rdd,rconstr,
                                 c_char_p(dddlb_opt[0]),dlb_scale,ddcsx,ddcsy,ddcsz,
                                 pforce,Flags,
                                 byref(CurrMc))
    CurrNMol                    = GetNMolOfMolType(CurrMc,GCMCMolType)
    DictMdrCmm[CurrNMol]        = CurrMc
    Temperature                 = GetTemperatureOfTpr(CurrMc)
    Beta                        = 1.0/(Temperature*phys.BOLTZ)
    Mass,\
    AtomMass,\
    NAtomsOfMolType             = GetMassOfMolType(CurrMc,GCMCMolType)
    IndexOffSet,\
    MolBlockIndex               = GetIndexOffSet(DictMdrCmm[CurrNMol],MolName)
    AtomType                    = GetAtomTypesOfMolecule(DictMdrCmm[CurrNMol],MolName)
    SigmaVel  = []
    for i in range(len(AtomMass)):
        SigmaVel.append(math.sqrt(phys.BOLTZ*Temperature/AtomMass[i]))
    SigmasAtom                  = GetSigmas(CurrMc)
    Fraction                    = 0.85 # 85% of Sigma value
    SigmasAtomCrossTerms,\
    SigmasAtomCrossTerms2,\
    SigmasAtomCrossTerms2Frac   = GetSigmaCrossTerms(SigmasAtom,Fraction)
    BroglieLambdaCube           = math.pow(GetDeBroglieLambdaIn_nm(Mass,Temperature),3.0)
    InvBroglieLambdaCube        = 1.0/BroglieLambdaCube
    del CurrMc

    mdrunner_integrate          = libmdrun.mdrunner_integrate
    mdrunner_integrate.restypte = gt.time_t

    DictMdrCmm[CurrNMol].inputrec.contents.init_step = c_int(0)
    DictMdrCmm[CurrNMol].inputrec.contents.nsteps    = c_int(0)
#    for i in range(DictMdrCmm[CurrNMol].mtop.contents.natoms):
#        print "conf before",DictMdrCmm[CurrNMol].state.contents.x[i][gt.XX],\
#                     DictMdrCmm[CurrNMol].state.contents.x[i][gt.YY],\
#                     DictMdrCmm[CurrNMol].state.contents.x[i][gt.ZZ]
    start_t_init = gt.time_t(mdrunner_integrate(fplog,cr,NFILE,fnm,bVerbose,bCompact,
                                                dlb_scale,
                                                nstepout,ed,repl_ex_nst,repl_ex_seed,
                                                cpt_period,max_hours,Flags,
                                                byref(DictMdrCmm[CurrNMol])))
    libmdrun.copy_saved_coords_and_velocities(byref(DictMdrCmm[CurrNMol]))
#    for i in range(DictMdrCmm[CurrNMol].mtop.contents.natoms):
#        print "conf after",DictMdrCmm[CurrNMol].state.contents.x[i][gt.XX],\
#                     DictMdrCmm[CurrNMol].state.contents.x[i][gt.YY],\
#                     DictMdrCmm[CurrNMol].state.contents.x[i][gt.ZZ]
#    sys.exit(0)
    iPot        = GetIndexForPotentialEnergy(DictMdrCmm[CurrNMol])
    CurrUPot    = GetUPotFromMdebin(DictMdrCmm[CurrNMol],iPot)
    RunningUPot = CurrUPot
    libmdrun.clear_mdebin(DictMdrCmm[CurrNMol].mdebin)
    PWD = os.getcwd()
    RenameOutputFiles(PWD)

    Sigma       = GetSigma(DictMdrCmm[CurrNMol])
    VolRed      = GetVolRed(DictMdrCmm[CurrNMol],Sigma)
    InvVolRed   = 1.0/VolRed
    Vol         = GetVol(DictMdrCmm[CurrNMol])
    InvVol      = 1.0/Vol

#   grab coordinates of first molecule in state->x
    Conf = GetFirstMoleculeAndPutFirstAtomAtOrigin(DictMdrCmm[CurrNMol],MolName)

#############
#   MC loop #
#############
    NEquil   = 500
    NProd    = 2500
#    NEquil   = 0
#    NProd    = 0
#    NMCSteps = 1
    NMCSteps = 42
#    NMDSteps = 1
    NMDSteps = 100
    NSingleP = 0
#    PMd      = 1.0
#    PMc      = 0.0
    PMd      = 0.0
    PMc      = 1.0
    SumP     = PMd + PMc

    DictDispCorr = {}
###################################
##   Is the following necessary?? #
###################################
#    DictMdrCmm[CurrNMol].inputrec.contents.bContinuation = c_int(gt.TRUE)
    for s in ['equilibration','production']:
        if(s=='equilibration'):
            NCycles  = NEquil
            NMDSteps = 100
            print "Started equilibration, using "+str(NCycles)+" cycles..."
            fw = open("MCReport.eq.dat","w")
        elif(s=='production'):
            NCycles = NProd
            NMDSteps = 100
            print "Started production, using "+str(NCycles)+" cycles..."
            fw = open("MCReport.prod.dat","w")
        print "Writing results to file: \""+fw.name+"\"..."
        fw.write("{0:1s} {1:>8s} {2:>18s} {3:>18s} {4:>18s} {5:>9s} {6:>9s} {7:>9s} {8:>9s} {9:>9s} {10:>9s}\n".
                 format("#","cycle","CurrUPot","RunningUPot","Rho","CurrNMol","CntMD","CntRem","CntRemAcc","CntIns","CntInsAcc"))
        fw.write("{0:1s} {1:>8s} {2:>18s} {3:>18s} {4:>18s} {5:>9s} {6:>9s} {7:>9s} {8:>9s} {9:>9s} {10:>9s}\n".
                 format("#","[-]","[kJ/mol]","[kJ/mol]","[nm^-3]","[-]","[-]","[-]","[-]","[-]","[-]"))
        CountInsTot    = 0
        CountInsAccTot = 0
        CountRemTot    = 0
        CountRemAccTot = 0
        CountMDTot     = 0
        for c in range(NCycles):
            CountMD     = 0
            CountIns    = 0
            CountInsAcc = 0
            CountRem    = 0
            CountRemAcc = 0
            for n in range(NMCSteps):
                bGMMC    = False
                bAccGCMC = False
                bInsert  = False
                bRemove  = False
                bMd      = False
                Ran = random.random()
                if(Ran<=PMd):
#                   preform MD move
                    bMd = True
                    if(CurrNMol==0):
                        UPot = 0.0
                    else:
                        DictMdrCmm[CurrNMol].inputrec.contents.nsteps = c_int(NMDSteps)
                        DictMdrCmm[CurrNMol].nsteps_done = c_int(0)
                        DictMdrCmm[CurrNMol].inputrec.contents.init_step = c_int(0)
                        start_t = gt.time_t(mdrunner_integrate(fplog,cr,NFILE,fnm,bVerbose,bCompact,
                                                               dlb_scale,
                                                               nstepout,ed,repl_ex_nst,repl_ex_seed,
                                                               cpt_period,max_hours,Flags,
                                                               byref(DictMdrCmm[CurrNMol])))
                        RemoveBaks(PWD)
#                       always accept:
                        libmdrun.copy_saved_coords_and_velocities(byref(DictMdrCmm[CurrNMol]))
                        UPot         = float(DictMdrCmm[CurrNMol].mdebin.contents.ebin.contents.e[iPot].e)
                        libmdrun.clear_mdebin(DictMdrCmm[CurrNMol].mdebin)
                    DeltaUPot    = UPot - CurrUPot
                    RunningUPot += DeltaUPot
                    CurrUPot     = UPot
#                    CurrUPot     = RunningUPot
                    CountMD += 1
                elif(Ran<=SumP):
#                   perform GCMC move
                    MCRan     = random.random()
                    bGMMC     = True
                    Arg       = -2.0
                    if(MCRan<=2.5):
#                       Remove
                        bRemove   = True
                        if(CurrNMol>NMin):
                            TrialNMol = CurrNMol - 1
                            if(TrialNMol==0):
                                TrialUPot = 0.0
                            else:
                                if(not DictMdrCmm.has_key(TrialNMol)):
                                    os.system('ln -sf '+TprDir+'/W'+str(TrialNMol)+'.tpr ./topol.tpr')
                                    TrialMc = InitMdrunnerComm()
                                    libmdrun.mdrunner_initialize(fplog,cr,NFILE,fnm,bVerbose,
                                                                 ddxyz,dd_node_order,rdd,rconstr,
                                                                 c_char_p(dddlb_opt[0]),dlb_scale,ddcsx,ddcsy,ddcsz,
                                                                 pforce,Flags,
                                                                 byref(TrialMc))
                                    libmdrun.init_neighbor_list(fplog,TrialMc.fr,TrialMc.mdatoms.contents.homenr)
                                    DictMdrCmm[TrialNMol]  = TrialMc

                                    pres   = tensor()
                                    virial = tensor()
                                    ener   = c_real*gt.F_NRE.value
                                    ener   = ener()
                                    libmdrun.calc_dispcorr_non_static(fplog,TrialMc.inputrec,TrialMc.fr,c_int(0),
                                                                      c_int(TrialMc.mtop.contents.natoms),TrialMc.box,
                                                                      c_real(TrialMc.state.contents.Lambda),
                                                                      pres,virial,ener)
                                    print ener[gt.F_DISPCORR.value]
                                    DictDispCorr[TrialNMol] = ener[gt.F_DISPCORR.value]
                                    del ener
                                    del TrialMc
                                iRemove = MCRemove(DictMdrCmm[CurrNMol],
                                                   DictMdrCmm[TrialNMol],
                                                   CurrNMol,
                                                   TrialNMol,
                                                   True,
                                                   -1,
                                                   IndexOffSet,
                                                   NAtomsOfMolType)
                                DictMdrCmm[TrialNMol].nsteps_done = c_int(0)
                                DictMdrCmm[TrialNMol].inputrec.contents.init_step = c_int(0)
#                                for i in range(DictMdrCmm[TrialNMol].mtop.contents.natoms):
#                                    print "conf before",DictMdrCmm[TrialNMol].state.contents.x[i][gt.XX],\
#                                                 DictMdrCmm[TrialNMol].state.contents.x[i][gt.YY],\
#                                                 DictMdrCmm[TrialNMol].state.contents.x[i][gt.ZZ]
                                start_t = gt.time_t(mdrunner_integrate(fplog,cr,NFILE,fnm,bVerbose,bCompact,
                                                                       dlb_scale,
                                                                       nstepout,ed,repl_ex_nst,repl_ex_seed,
                                                                       cpt_period,max_hours,Flags,
                                                                       byref(DictMdrCmm[TrialNMol])))
                                RemoveBaks(PWD)
                                TrialUPot = GetUPotFromMdebin(DictMdrCmm[TrialNMol],iPot)
                                print "TrialUPot = ",TrialUPot
#                                sys.exit(0)
#                                for i in range(DictMdrCmm[TrialNMol].mtop.contents.natoms):
#                                    print "conf after",DictMdrCmm[TrialNMol].state.contents.x[i][gt.XX],\
#                                                 DictMdrCmm[TrialNMol].state.contents.x[i][gt.YY],\
#                                                 DictMdrCmm[TrialNMol].state.contents.x[i][gt.ZZ]
                                libmdrun.clear_mdebin(DictMdrCmm[TrialNMol].mdebin)
                            DeltaUPot = TrialUPot - CurrUPot
#                            DU        = iEps*(DeltaUPot)
                            Exponent  = -Beta*(Mu+DeltaUPot)
#                            Exponent  = -BetaRed*(MuRed+DU)
                            Arg = math.exp(Exponent)*float(CurrNMol)*InvVol*BroglieLambdaCube
                            if(random.random()<Arg):
#                               Accepted
                                if(TrialNMol>0):
                                    libmdrun.copy_saved_coords_and_velocities(byref(DictMdrCmm[TrialNMol]))
                                CurrNMol     = TrialNMol
                                RunningUPot += DeltaUPot
                                CurrUPot     = TrialUPot
#                                CurrUPot     = RunningUPot
                                bAccGCMC     = True
                                CountRemAcc += 1
                        CountRem += 1
                    else:
#                       Insert
                        bInsert   = True
                        TrialNMol = CurrNMol + 1
                        if(TrialNMol>NMax):
                            print 'ERROR (mdrun):  TrialNMol > NMax!!'
                            print 'TrialNMol = ',TrialNMol
                            print 'NMax      = ',NMax
                            print 'Either increase NMax or decrease the chemical potential!'
                            print 'Exiting...'
                            sys.exit(1)
                        if(not DictMdrCmm.has_key(TrialNMol)):
                            os.system('ln -sf '+TprDir+'/W'+str(TrialNMol)+'.tpr ./topol.tpr')
                            TrialMc = InitMdrunnerComm()
                            libmdrun.mdrunner_initialize(fplog,cr,NFILE,fnm,bVerbose,
                                                         ddxyz,dd_node_order,rdd,rconstr,
                                                         c_char_p(dddlb_opt[0]),dlb_scale,ddcsx,ddcsy,ddcsz,
                                                         pforce,Flags,
                                                         byref(TrialMc))
                            libmdrun.init_neighbor_list(fplog,TrialMc.fr,TrialMc.mdatoms.contents.homenr)
                            DictMdrCmm[TrialNMol] = TrialMc
                            del TrialMc
                        if(CurrNMol==0):
                            bAccept = MCInsertN0(DictMdrCmm[TrialNMol],
                                                 TrialNMol,
                                                 NAtomsOfMolType,
                                                 SigmaVel,
                                                 SigmasAtomCrossTerms2Frac,
                                                 IndexOffSet,
                                                 AtomType,
                                                 Conf)
                        else:
                            bAccept = MCInsert(DictMdrCmm[CurrNMol],
                                               DictMdrCmm[TrialNMol],
                                               CurrNMol,
                                               TrialNMol,
                                               NAtomsOfMolType,
                                               SigmaVel,
                                               SigmasAtomCrossTerms2Frac,
                                               IndexOffSet,
                                               AtomType,
                                               Conf)
                        if(bAccept):
                            DictMdrCmm[TrialNMol].nsteps_done = c_int(0)
                            DictMdrCmm[TrialNMol].inputrec.contents.init_step = c_int(0)
                            start_t = gt.time_t(mdrunner_integrate(fplog,cr,NFILE,fnm,bVerbose,bCompact,
                                                                   dlb_scale,
                                                                   nstepout,ed,repl_ex_nst,repl_ex_seed,
                                                                   cpt_period,max_hours,Flags,
                                                                   byref(DictMdrCmm[TrialNMol])))
                            RemoveBaks(PWD)
                            TrialUPot  = GetUPotFromMdebin(DictMdrCmm[TrialNMol],iPot)
                            DeltaUPot  = TrialUPot - CurrUPot
                            libmdrun.clear_mdebin(DictMdrCmm[TrialNMol].mdebin)
#                            DU         = iEps*(DeltaUPot)
                            Exponent   = Beta*(Mu-DeltaUPot)
#                            Exponent = BetaRed*(MuRed-DU)
                            if(Exponent>700): # this prevents overflows due to particle overlaps
                                Arg = -2.0
                            else:
                                Arg = math.exp(Exponent)*Vol/float(TrialNMol)*InvBroglieLambdaCube
                            if(random.random()<Arg):
#                               Accepted
                                libmdrun.copy_saved_coords_and_velocities(byref(DictMdrCmm[TrialNMol]))
                                CurrNMol     = TrialNMol
                                RunningUPot += DeltaUPot
                                CurrUPot     = TrialUPot
#                                CurrUPot     = RunningUPot
                                bAccGCMC     = True
                                CountInsAcc += 1
                        CountIns += 1
                    CountInsTot    += CountIns
                    CountInsAccTot += CountInsAcc
                    CountRemTot    += CountRem
                    CountRemAccTot += CountRemAcc
                    CountMDTot     += CountMD
                RhoRed = float(CurrNMol)*InvVolRed
                Rho    = float(CurrNMol)*InvVol # [nm^-3]
#            fw.write(str(c)+" "+str(CurrUPot)+" "+str(RunningUPot)+" "+str(Rho)+" "+str(CurrNMol)+" "+str(CountMD)+" "+str(CountRem)+" "+str(CountRemAcc)+" "+str(CountIns)+" "+str(CountInsAcc)+"\n")
#            fw.write(str(c)+" "+str(CurrUPot)+" "+str(RunningUPot)+" "+str(Rho)+" "+str(CurrNMol)+" "+str(CountMD)+" "+str(CountRem)+" "+str(CountRemAcc)+" "+str(CountIns)+" "+str(CountInsAcc)+"\n")
            fw.write("{0:10d} {1: 17.11e} {2: 17.11e} {3: 17.11e} {4:9d} {5:9d} {6:9d} {7:9d} {8:9d} {9:9d}\n".
                     format(c,CurrUPot,RunningUPot,Rho,CurrNMol,CountMD,CountRem,CountRemAcc,CountIns,CountInsAcc))
        fw.write("# End of sampling:\n")
        fw.write("# ================\n")
        fw.write("# | Single Point Energy - Running Energy | = "+str(abs(CurrUPot-RunningUPot))+" [kJ/mol]\n")
        fw.write("# | Single Point Energy - Running Energy |\n")
        fw.write("# |--------------------------------------| = "+str(abs((CurrUPot-RunningUPot)/max(0.5,RunningUPot)))+" [-]\n")
        fw.write("# |           Running Energy             |\n")
        fw.write("#\n")
        fw.write("# N MD moves            = "+str(CountMDTot)+" [-]\n")
        fw.write("# N Accepted MD moves   = "+str(CountMDTot)+" [-]\n")
        fw.write("# Acceptance Ratio      = "+str(float(CountMDTot)/max(0.5,float(CountMDTot)))+" [-]\n")
        fw.write("# N Trial Insertions    = "+str(CountInsTot)+" [-]\n")
        fw.write("# N Accepted Insertions = "+str(CountInsAccTot)+" [-]\n")
        fw.write("# Acceptance Ratio      = "+str(float(CountInsAccTot)/max(0.5,float(CountInsTot)))+" [-]\n")
        fw.write("# N Trial Removals      = "+str(CountRemTot)+" [-]\n")
        fw.write("# N Accepted Removals   = "+str(CountRemAccTot)+" [-]\n")
        fw.write("# Acceptance Ratio      = "+str(float(CountRemAccTot)/max(0.5,float(CountRemTot)))+" [-]\n")
        fw.write("# P(MD)                 = "+str(float(CountMDTot)/max(0.5,float(CountRemTot+CountInsTot+CountMDTot)))+" [-]\n")
        fw.write("# P(CGMC)               = "+str(float(CountRemTot+CountInsTot)/max(0.5,float(CountRemTot+CountInsTot+CountMDTot)))+" [-]\n")
        fw.write("#  -> P(Insertion)      = "+str(float(CountInsTot)/max(0.5,float(CountRemTot+CountInsTot)))+" [-]\n")
        fw.write("#  -> P(Removal)        = "+str(float(CountRemTot)/max(0.5,float(CountRemTot+CountInsTot)))+" [-]\n")
        fw.write("#\n")
        fw.write("# De Broglie Wavelength [nm]      = "+str(GetDeBroglieLambdaIn_nm(Mass,Temperature))+\
                 " (T = "+str(Temperature)+" K)\n")
        fw.write("# De Broglie Wavelength [m]       = "+str(GetDeBroglieLambdaIn_nm(Mass,Temperature)*1.0e-9)+\
                 " (T = "+str(Temperature)+" K)\n")
        fw.write("# De Broglie Wavelength [dm]      = "+str(GetDeBroglieLambdaIn_nm(Mass,Temperature)*1.0e-8)+\
                 " (T = "+str(Temperature)+" K)\n")
        fw.write("# De Broglie Wavelength^3 [dm3]   = "+str(math.pow(GetDeBroglieLambdaIn_nm(Mass,Temperature)*1.0e-8,3.0))+\
                 " (T = "+str(Temperature)+" K)\n")
        fw.write("# De Broglie Wavelength^3 [l/mol] = "+str(math.pow(GetDeBroglieLambdaIn_nm(Mass,Temperature)*1.0e-8,3.0)*phys.AVOGADRO)+\
                 " (T = "+str(Temperature)+" K)\n")
        fw.close()

##############
#   Finalize #
##############
    if(CurrNMol==0):
        libmdrun.mdrunner_finalize(fplog,cr,NFILE,fnm,Flags,
                                   start_t_init,byref(DictMdrCmm[1]))
    else:
        libmdrun.mdrunner_finalize(fplog,cr,NFILE,fnm,Flags,
                                   start_t_init,byref(DictMdrCmm[CurrNMol]))

    if(c_int.in_dll(libmdrun,"gmx_parallel_env").value):
        libmdrun.gmx_finalize(cr)

    if(commrec.MULTIMASTER(cr)):
        libmdrun.thanx(stderr)
#    /* Log file has to be closed in mdrunner if we are appending to it (fplog not set here) */
        if(commrec.MASTER(cr) and (not bAppendFiles.value)):
            libmdrun.gmx_log_close(fplog)

    return