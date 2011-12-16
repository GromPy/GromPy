#! /bin/bash

for GroFile in `ls ../gro | grep \.gro`
do
	BaseName=`echo $GroFile | cut -d. -f1`
	grompp -f ../mdp/gcmc.mdp -c ../gro/$GroFile -p ../top/$BaseName.top -o $BaseName.tpr -maxwarn 1
	rm -f \#*
    rm mdout.mdp
done

exit
