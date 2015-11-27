#!/bin/sh

if [ -z $1 ]; then
    echo "Usage: $0 [bool electron] [bool muon]"
    exit 0
fi

electron=$1
muon=$2

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src
cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

if [ $electron == true ]; then

    root -q -b -l dataMCplots.C++\(\"electron/dataMC_Ele/outputEle/barrel\"\,\"barrelEleVariables\"\)
    root -q -b -l dataMCplots.C++\(\"electron/dataMC_Ele/outputEle/endcap\"\,\"endcapEleVariables\"\)
    root -q -b -l dataMCplots.C++\(\"electron/dataMC_Zee/outputZee\"\,\"ZeeVariables\"\)
    root -q -b -l dataMCplots.C++\(\"electron/dataMC_Jetee/outputJetee\"\,\"JeteeVariables\"\)

else

    electron=false

fi

if [ $muon == true ]; then

    root -q -b -l dataMCplots.C++\(\"muon/dataMC_Mu/outputMu/highpt\"\,\"highPtMuVariables\"\)
    root -q -b -l dataMCplots.C++\(\"muon/dataMC_Mu/outputMu/tracker\"\,\"trackerMuVariables\"\)
    root -q -b -l dataMCplots.C++\(\"muon/dataMC_Zmumu/outputZmumu\"\,\"ZmumuVariables\"\)
    root -q -b -l dataMCplots.C++\(\"muon/dataMC_Jetmumu/outputJetmumu\"\,\"JetmumuVariables\"\)
    root -q -b -l dataMCplots.C++\(\"muon/forLimitCalc/outputmZHmu\"\,\"mZHmuVariables\"\)

else

    muon=false

fi

rm -f *.d *.so
mv *pdf /afs/cern.ch/user/h/htong/www

exit