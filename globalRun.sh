#!/bin/sh

if [ -z $6 ]; then

    echo "Usage: $0 [macro without .C] [bool data] [bool DY] [bool diboson] [bool ttbar] [bool signal]"
    exit 0

fi

data=$2
DY=$3
diboson=$4
ttbar=$5
signal=$6

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src
cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

mudatapath=/data7/syu/NCUGlobalTuples/Run2015D/9b33d00/SingleMuon
eledatapath=/data7/khurana/NCUGlobalTuples/Run2015DFullDataset

if [ $data == true ]; then

    if [ `echo $1 | grep -c "ele"` -gt 0 ]; then

	root -q -b -l $1.C++\(\"$eledatapath/crab_SingleElectron-Run2015D-05Oct2015-v1_20151117_2p2fb_SingleEleTextFile/\"\,\"SingleElectron-Run2015D-V120151117\"\)
	root -q -b -l $1.C++\(\"$eledatapath/crab_SingleElectron-Run2015D-PromptReco-V420151117_2p2fb_SingleEleTextFile/\"\,\"SingleElectron-Run2015D-V420151117\"\)

    elif [ `echo $1 | grep -c "mu"` -gt 0 ]; then

	root -q -b -l $1.C++\(\"$mudatapath/crab_SingleMuon-Run2015D-05Oct2015-v1_20151119_2p2fb_SingleMuTextFile/\"\,\"SingleMuon-Run2015D-V120151119\"\)             
	root -q -b -l $1.C++\(\"$mudatapath/crab_SingleMuon-Run2015D-PromptReco-V420151119_2p2fb_SingleMuTextFile/\"\,\"SingleMuon-Run2015D-V420151119\"\)

    fi

else

    data=false

fi

mcpath=/data7/syu/NCUGlobalTuples/Spring15_ReMiniAODSim/9b33d00

if [ $DY == true ]; then

    root -q -b -l $1.C++\(\"$mcpath/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\)
    root -q -b -l $1.C++\(\"$mcpath/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\)
    root -q -b -l $1.C++\(\"$mcpath/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\)
    root -q -b -l $1.C++\(\"$mcpath/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\)

else

    DY=false

fi

if [ $diboson == true ]; then

    root -q -b -l $1.C++\(\"$mcpath/WW_TuneCUETP8M1_13TeV-pythia8/\"\,\"WW_TuneCUETP8M1_13TeV\"\)
    root -q -b -l $1.C++\(\"$mcpath/WZ_TuneCUETP8M1_13TeV-pythia8/\"\,\"WZ_TuneCUETP8M1_13TeV\"\)
    root -q -b -l $1.C++\(\"$mcpath/ZZ_TuneCUETP8M1_13TeV-pythia8/\"\,\"ZZ_TuneCUETP8M1_13TeV\"\)

else

    diboson=false

fi

if [ $ttbar == true ]; then

    root -q -b -l $1.C++\(\"$mcpath/TT_TuneCUETP8M1_13TeV-powheg-pythia8/\"\,\"TT_TuneCUETP8M1_13TeV\"\)

else

    ttbar=false

fi

if [ $signal == true ]; then

    mass=(600 800 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000 4500)

    for i in "${mass[@]}"
    do

	sigalname=ZprimeToZhToZlephbb_narrow_M-${i}_13TeV-madgraph
	root -q -b -l $1.C++\(\"$mcpath/ZprimeToZhToZlephbb/$sigalname/\"\,\"ZprimeToZhToZlephbb_M-${i}_13TeV\"\)

    done
   
else

    signal=false

fi

rm -f inputdir.txt
rm -f *.pcm *.d *.so
mv *root output*

exit