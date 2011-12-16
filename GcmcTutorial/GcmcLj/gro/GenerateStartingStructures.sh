#! /bin/bash

EXPECTED_N_ARGS=3
if [ $# -ne $EXPECTED_N_ARGS ]; then
   echo "Input the correct amount of arguments!"
   echo "Usage:"
   echo "------"
   echo "./GenerateStartingStructures.sh \"StartOfNRange\" \"EndOfNRange\" \"WInputN.gro\""
   echo "Example: ./GenerateStartingStructures.sh 1 430 W400.gro"
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
    head -1 $NInputFile                     > tmp.gro
    echo $N                                >> tmp.gro
    head -$(($N+2)) $NInputFile | tail -$N >> tmp.gro
    tail -1 $NInputFile                    >> tmp.gro
    mv tmp.gro W$N.gro
    N=$N+1
done
declare -i N=$(($NInput+1))
while [ $N -le $NEnd ]
do
    PrevN=$(($N-1))
    head -1 $NInputFile                                        > tmp.gro
    echo $N                                                   >> tmp.gro
    head -$(($PrevN+2)) W$PrevN.gro | tail -$PrevN            >> tmp.gro
    tail -2 W$PrevN.gro | head -1 | sed -e s/"$PrevN"W/"$N"W/ >> tmp.gro
    tail -1 $NInputFile                                       >> tmp.gro
    editconf -f tmp.gro -o W$N.gro
    rm -f tmp.gro
    rm -f \#*
    N=$N+1
done

exit
