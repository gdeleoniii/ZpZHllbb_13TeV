#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

channel=(ele mu)
cat=(1 2)

for ((i=0; i<${#channel[@]}; i++)); do
    for ((j=0; j<${#cat[@]}; j++)); do
	echo "Now running ==> channel: " ${channel[$i]} " category: " ${cat[$j]}
	root -q -b -l rooFitData.C\(\"${channel[$i]}\"\,\"${cat[$j]}\"\,false\)
    done
done

#mkdir dataJetEnScaleResults/
#mv *pdf dataJetEnScaleResults/
#rm -rf $HOME/www/dataJetEnScaleResults/
#mv dataJetEnScaleResults/ $HOME/www/

echo "All the jobs are finished."

exit
