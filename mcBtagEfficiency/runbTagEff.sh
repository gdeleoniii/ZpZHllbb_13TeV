#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

channel=(ele mu)
flavor=(udsg c b)

for ((i=0; i<${#channel[@]}; i++)); do
    for ((k=0; k<${#flavor[@]}; k++)); do
	    root -q -b -l bTagEff.C+\(\"${channel[$i]}\"\,\"${flavor[$k]}\"\)
    done
done

mkdir MCbTagEfficiency/
mv *pdf MCbTagEfficiency/
rm -rf $HOME/www/MCbTagEfficiency/
mv MCbTagEfficiency/ $HOME/www/

echo "All the jobs are finished."

exit
