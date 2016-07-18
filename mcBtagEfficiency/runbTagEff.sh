#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

channel=(ele mu)
cat=(1 2)
flavor=(1 4 5)

for ((i=0; i<${#channel[@]}; i++)); do
    for ((j=0; j<${#cat[@]}; j++)); do
#	for ((k=0; k<${#flavor[@]}; k++)); do
#	    echo "Now running ==> channel: " ${channel[$i]} " category: " ${cat[$j]} " flavor: " ${flavor[$k]}
	    root -q -b -l bTagEff.C\(\"${channel[$i]}\"\,${cat[$j]}\,0\)
#	done
    done
done

mkdir MCbTagEfficiency/
mv *pdf MCbTagEfficiency/
rm -rf $HOME/www/MCbTagEfficiency/
mv MCbTagEfficiency/ $HOME/www/

echo "All the jobs are finished."

exit
