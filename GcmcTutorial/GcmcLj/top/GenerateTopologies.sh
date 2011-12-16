#! /bin/bash

EXPECTED_N_ARGS=3
if [ $# -ne $EXPECTED_N_ARGS ]; then
   echo "Input the correct amount of arguments!"
   echo "Usage:"
   echo "------"
   echo "./GenerateTopologies.sh \"StartOfNRange\" \"EndOfNRange\" \"WInputN.top\""
   echo "Example: ./GenerateTopologies.sh1 430 W400.top"
   echo "EXITING..."
   exit 1
fi

NStart=$1
NEnd=$2
NInputFile=$3
NInput=`echo $NInputFile | cut -d. -f1 | sed -e s/W//`

declare -i N=$NStart
while [ $N -lt $NInput ]
do
    NLines=`cat $NInputFile | wc -l`
    head -$(($NLines-1)) $NInputFile            > tmp.top
    tail -1 $NInputFile | sed -e s/$NInput/$N/ >> tmp.top
    mv tmp.top W$N.top
    N=$N+1
done
declare -i N=$(($NInput+1))
while [ $N -le $NEnd ]
do
    NLines=`cat $NInputFile | wc -l`
    head -$(($NLines-1)) $NInputFile            > tmp.top
    tail -1 $NInputFile | sed -e s/$NInput/$N/ >> tmp.top
    mv tmp.top W$N.top
    N=$N+1
done

exit
