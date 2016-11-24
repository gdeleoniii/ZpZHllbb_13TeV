#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

pwd=$PWD
#processEle=(Barrel Endcap Zjet Higgs)
#processMu=(HighPt Tracker Zjet Higgs)

processEle=Higgs
processMu=Higgs

cd $pwd/ele
for ((j=0; j<${#processEle[@]}; j++)); do

#    /bin/bash $pwd/globalRun.sh ele${processEle[$j]}Variable ele 1 1 1 1 1 0
#    mkdir -p output_ele${processEle[$j]}Variable
#    mv *root output_ele${processEle[$j]}Variable
    root -q -b -l $pwd/dataMCplots.C+\(\"Electron\"\,\"output_ele${processEle[$j]}Variable\"\,\"ele${processEle[$j]}Variable\"\)

done

rm -f *.pcm *.d *.so
mv *pdf /afs/cern.ch/user/h/htong/www/dataMCcompare

cd $pwd/mu
for ((j=0; j<${#processMu[@]}; j++)); do

#    /bin/bash $pwd/globalRun.sh mu${processMu[$j]}Variable mu 1 1 1 1 1 0
#    mkdir -p output_mu${processMu[$j]}Variable
#    mv *root output_mu${processMu[$j]}Variable
    root -q -b -l $pwd/dataMCplots.C+\(\"Muon\"\,\"output_mu${processMu[$j]}Variable\"\,\"mu${processMu[$j]}Variable\"\)

done

rm -f *.pcm *.d *.so
mv *pdf /afs/cern.ch/user/h/htong/www/dataMCcompare

cd $pwd
rm -f *.pcm *.d *.so

exit