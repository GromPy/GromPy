#! /bin/bash

MCReportFiles=`find -type f -name MCReport.prod.dat`

echo "# mu [kJ/mol] Rho [1/(nm^3)] Err Rho [1/(nm^3)]" > rho_of_mu.dat
for MCReportFile in $MCReportFiles
do
    grep -v \# $MCReportFile | awk '{print $1,$4}' > Rho.xvg
    g_analyze -f Rho.xvg -ee Err.xvg > /dev/null 2>&1
    rm -f \#*
    Mu=`echo $MCReportFile | sed -e s/"\/"/\ /g | awk '{print $2}' | sed -e s/mu//`
    AvRho=`grep \"av Err.xvg | sed -e s/\"//g | awk '{print $(NF)}'`
    SdRho=`grep \"ee Err.xvg | sed -e s/\"//g | awk '{print $(NF)}'`
    echo $Mu $AvRho $SdRho
done | sort -n >> rho_of_mu.dat
rm Rho.xvg Err.xvg

NA="6.0221415e23"
MCReportFile=`echo $MCReportFiles | awk '{print $(NF)}'`
Lambda3=`grep "\# De Broglie Wavelength\^3 \[l\/mol\]" $MCReportFile | awk '{print $(NF-4)}'` #l/mol
DmPerM=10.0
NmPerM=1.0e9
Kilo=1000.0
R=8.314472 #J/mol/K
T=`grep "\# De Broglie Wavelength\^3 \[l\/mol\]" $MCReportFile | awk '{print $(NF-1)}'` #K

echo "#  mu [kJ/mol] Rho [mol/l] Err Rho [mol/l]"                                   > rho_of_mu_in_mol_p_l.dat
grep -v \# rho_of_mu.dat | awk '{f='$NmPerM'^3/'$Kilo'/'$NA';print $1,f*$2,f*$3}'  >> rho_of_mu_in_mol_p_l.dat

echo "# rho [mol/l] sdrho [mol/l] mu [kJ/mol] muid[kJ/mol] muex [kJ/mol] sdmuex [kJ/mol]" > mu_of_rho_in_mol_p_l.dat
grep -v \# rho_of_mu_in_mol_p_l.dat | awk '
{
rho=$2;
sdrho=$3;
mu=$1;
muid='$R'*'$T'*log(rho*'$Lambda3')/'$Kilo';
muex=mu-muid;
sdmuex='$R'*'$T'/rho*sdrho/'$Kilo';
print rho, sdrho, mu, muid, muex, sdmuex;
}' >> mu_of_rho_in_mol_p_l.dat

exit
