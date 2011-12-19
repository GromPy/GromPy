#! /bin/bash

GMXINSTALLDIR="../../gromacs/gromacs-4.0.7/gromacs-4.0.7-git/install/"
cd $GMXINSTALLDIR
GMXINSTALLDIR=$PWD
cd -
source ../SourceGromPyEnv.sh

InputParFile="start.dat"
MuList=`grep -v \# $InputParFile | awk '{print $3}' | sort -n`
TprDir=$PWD"/tpr"
PyLogFile="py.log"
NMin=0
NMax=426
GCMCMol="W"
GcmcPy=testHybrid.py

for i in $MuList
do
    Mu=$i
    NStart=`grep "mu"$i $InputParFile| awk '{print $2}'`
    echo "Performing MuVT simulation @ Mu=$Mu ..."
    cd "mu$i"
        time ../../../$GcmcPy $TprDir $Mu $NStart $GCMCMol $NMin $NMax 2> /dev/null > $PyLogFile
        rm -f \#*
        gzip -f md.log
        gzip -f py.log
        echo
    cd ../
done

exit
