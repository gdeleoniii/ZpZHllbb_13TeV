#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

channel=(ele mu)

for ((i=0; i<${#channel[@]}; i++)); do
    root -q -b -l bTaggingEff.C+\(\"${channel[$i]}\"\)
done

rm -f *.pcm *.d *.so

mkdir signalBtagEfficiency/
mv *pdf signalBtagEfficiency/
rm -rf $HOME/www/signalBtagEfficiency/
mv signalBtagEfficiency/ $HOME/www/

echo "All the jobs are finished."

exit
