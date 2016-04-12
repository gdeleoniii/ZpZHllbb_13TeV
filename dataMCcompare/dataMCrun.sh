#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src
cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

processEle=(Barrel Endcap Z Jet)
processMu=(HighPt Tracker Z Jet)

cd $pwd/ele

for ((j=0; j<${#processEle[@]}; j++)); do

    thisName=ele${processEle[$j]}Variable

    /bin/bash $pwd/globalRun.sh $thisName ele true true true true false false

    mv *root output_$thisName

    root -q -b -l $pwd/dataMCplots.C+\(\"output_$thisName\"\,\"$thisName\"\)

done

rm -f *.pcm *.d *.so
mv *pdf /afs/cern.ch/user/h/htong/www/dataMCcompare_v2

cd $pwd/mu

for ((j=0; j<${#processMu[@]}; j++)); do

    thisName=mu${processMu[$j]}Variable

    /bin/bash $pwd/globalRun.sh $thisName mu true true true true false false

    mv *root output_$thisName

    root -q -b -l $pwd/dataMCplots.C+\(\"output_$thisName\"\,\"$thisName\"\)

done

rm -f *.pcm *.d *.so
mv *pdf /afs/cern.ch/user/h/htong/www/dataMCcompare_v2

cd $pwd
rm -f *.pcm *.d *.so

exit