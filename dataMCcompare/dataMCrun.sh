#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_4_15/src
cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

processEle=(Barrel Endcap Z Jet)
cd $pwd/ele

for ((j=0; j<${#processEle[@]}; j++)); do

#    /bin/bash $pwd/globalRun.sh ele${processEle[$j]}Variable ele 1 1 1 1 1 0
#    mv *root output_ele${processEle[$j]}Variable
    root -q -b -l $pwd/dataMCplots.C+\(\"Electron\"\,\"output_ele${processEle[$j]}Variable\"\,\"ele${processEle[$j]}Variable\"\)

done

rm -f *.pcm *.d *.so
mv *pdf /afs/cern.ch/user/h/htong/www/dataMCcompare

processMu=(HighPt Tracker Z Jet)
cd $pwd/mu

for ((j=0; j<${#processMu[@]}; j++)); do

#    /bin/bash $pwd/globalRun.sh mu${processMu[$j]}Variable mu 1 1 1 1 1 0
#    mv *root output_mu${processMu[$j]}Variable
    root -q -b -l $pwd/dataMCplots.C+\(\"Muon\"\,\"output_mu${processMu[$j]}Variable\"\,\"mu${processMu[$j]}Variable\"\)

done

rm -f *.pcm *.d *.so
mv *pdf /afs/cern.ch/user/h/htong/www/dataMCcompare

cd $pwd
rm -f *.pcm *.d *.so

exit