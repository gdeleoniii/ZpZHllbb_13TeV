#!/bin/sh

if [ -z $7 ]; then

    echo "Usage: $0 [macro without .C] [channel] [data] [DY] [diboson] [ttbar] [signal]"
    exit 0

fi

channel=$2
data=$3
DY=$4
diboson=$5
ttbar=$6
signal=$7

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src
cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

samplePath=/data7/htong/skim_samples/$channel

if [ $data == true ]; then

    if [ $channel == ele ]; then

	root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_SingleElectron-Run2015D-05Oct2015-v1_20151117_2p2fb_SingleEleTextFile.root\"\,\"SingleElectron-Run2015D-V120151117\"\)
	root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_SingleElectron-Run2015D-PromptReco-V420151117_2p2fb_SingleEleTextFile.root\"\,\"SingleElectron-Run2015D-V420151117\"\)

    elif [ $channel == mu ]; then

	root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_SingleMuon-Run2015D-05Oct2015-v1_20151119_2p2fb_SingleMuTextFile.root\"\,\"SingleMuon-Run2015D-V120151119\"\)             
	root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_SingleMuon-Run2015D-PromptReco-V420151119_2p2fb_SingleMuTextFile.root\"\,\"SingleMuon-Run2015D-V420151119\"\)

    fi

else

    data=false

fi

if [ $DY == true ]; then

    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\)
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\)
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\)
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\)

else

    DY=false

fi

if [ $diboson == true ]; then

    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_WW_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WW_TuneCUETP8M1_13TeV\"\)
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_WZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WZ_TuneCUETP8M1_13TeV\"\)
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_ZZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"ZZ_TuneCUETP8M1_13TeV\"\)

else

    diboson=false

fi

if [ $ttbar == true ]; then

    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root\"\,\"TT_TuneCUETP8M1_13TeV\"\)

else

    ttbar=false

fi

if [ $signal == true ]; then

    mass=(600 800 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000 4500)

    for i in "${mass[@]}"
    do

	signalName=skim_${channel}_ZprimeToZhToZlephbb_narrow_M-${i}_13TeV-madgraph.root
	root -q -b -l $1.C+\(\"$samplePath/$signalName\"\,\"ZprimeToZhToZlephbb_M-${i}_13TeV\"\)

    done
   
else 

    signal=false

fi

rm -f inputdir.txt
rm -f *.pcm *.d *.so

exit