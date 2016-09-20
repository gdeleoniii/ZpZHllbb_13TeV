#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

channel=(ele mu)
cat=(1 2)

for ((i=0; i<${#channel[@]}; i++)); do
    for ((j=0; j<${#cat[@]}; j++)); do
	echo "Now running ==> channel: " ${channel[$i]} " category: " ${cat[$j]}
	root -q -b -l getAlphaUnc.C\(\"${channel[$i]}\"\,\"${cat[$j]}\"\,\"mur1\"\,0\,2\)
	root -q -b -l getAlphaUnc.C\(\"${channel[$i]}\"\,\"${cat[$j]}\"\,\"pdf\"\,9\,109\)
    done
done

mkdir alphaPdfScaleResults/
mv *pdf alphaPdfScaleResults/
rm -rf $HOME/www/alphaPdfScaleResults/
mv alphaPdfScaleResults/ $HOME/www/

echo "All the jobs are finished."

exit
