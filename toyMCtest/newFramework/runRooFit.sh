#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

channel=(ele mu)
phy=(Zjets VV TT)
cat=(1 2)

for ((i=0; i<${#channel[@]}; i++)); do
    for ((j=0; j<${#cat[@]}; j++)); do
	for ((k=0; k<${#phy[@]}; k++)); do
	    echo "Now running ==> channel: " ${channel[$i]} " category: " ${cat[$j]} " physics process: " ${phy[$k]}
	    root -q -b -l rooFitTest.C\(\"${channel[$i]}\"\,\"${phy[$k]}\"\,\"${cat[$j]}\"\)
	done
    done
done

mkdir rooFittoyMCResults/
mv *pdf rooFittoyMCResults/

echo "All the jobs are finished."

exit
