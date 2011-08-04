#! /bin/bash

InputParFile="start.dat"
MuList=`grep -v \# $InputParFile | awk '{print $3}' | sort -n`
TprDir=$PWD"/tpr"
PyLogFile="py.log"
NMin=0
NMax=426
GCMCMol="W"

for i in $MuList
do
  Mu=$i
  NStart=`grep "mu"$i $InputParFile| awk '{print $2}'`
  echo "Performing MuVT simulation @ Mu=$Mu ..."
  cd "mu$i"
  time python /home/rpool/workspace/PyWork/trunk/GromPy/testHybrid.py $TprDir $Mu $NStart $GCMCMol $NMin $NMax 2> /dev/null > $PyLogFile
  gzip -f md.log
  gzip -f py.log 
  echo
  cd ../
done

exit
