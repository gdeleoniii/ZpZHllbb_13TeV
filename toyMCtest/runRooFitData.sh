#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

channel=(ele mu)

for ((i=0; i<${#channel[@]}; i++)); do
    echo "Now running ==> channel: " ${channel[$i]}
    root -q -b -l rooFitData.C\(\"${channel[$i]}\"\)
done

mkdir rooFitDataResults/
mv *pdf rooFitDataResults/
rm -rf $HOME/www/rooFitDataResults/
mv rooFitDataResults/ $HOME/www/

echo "All the jobs are finished."

exit
