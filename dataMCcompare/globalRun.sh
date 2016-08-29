#!/bin/sh

if [ -z $8 ]; then

    echo "Usage: $0 [macro without .C] [channel] [data] [DY] [diboson] [zh] [ttbar] [signal]"
    exit 0

fi

channel=$2
data=$3
DY=$4
diboson=$5
zh=$6
ttbar=$7
signal=$8
pwd=$PWD

#Do cmsenv:
#cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src
#cd $cmsswdr
#export SCRAM_ARCH=slc6_amd64_gcc481
#eval `scramv1 runtime -sh`
#cd $pwd

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

samplePath=/data7/htong/skim_NCUGlobalTuples

if [ $data == 1 ]; then

    if [ $channel == ele ]; then

	root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_SingleElectron_Run2015D-05Oct2015-v1.root\"\,\"SingleElectron-Run2015D-v1\"\)
	root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_SingleElectron_Run2015D-PromptReco-v4.root\"\,\"SingleElectron-Run2015D-v4\"\)

    elif [ $channel == mu ]; then

	root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_SingleMuon_Run2015D-05Oct2015-v1.root\"\,\"SingleMuon-Run2015D-v1\"\)             
	root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_SingleMuon_Run2015D-PromptReco-v4.root\"\,\"SingleMuon-Run2015D-v4\"\)

    fi

fi

if [ $DY == 1 ]; then
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\)
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\)
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\)
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\)
fi

if [ $diboson == 1 ]; then
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_WW_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WW_TuneCUETP8M1_13TeV\"\)
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_WZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WZ_TuneCUETP8M1_13TeV\"\)
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_ZZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"ZZ_TuneCUETP8M1_13TeV\"\)
fi

if [ $zh == 1 ]; then
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8.root\"\,\"ZH_HToBB_ZToLL_M125_13TeV\"\)
fi

if [ $ttbar == 1 ]; then
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root\"\,\"TT_TuneCUETP8M1_13TeV\"\)
fi

if [ $signal == 1 ]; then

    mass=(600 800 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000 4500)

    for i in "${mass[@]}"
    do

	signalName=skim_${channel}_crab_ZprimeToZhToZlephbb_narrow_M-${i}_13TeV-madgraph.root
	root -q -b -l $1.C+\(\"$samplePath/$signalName\"\,\"ZprimeToZhToZlephbb_M-${i}_13TeV\"\)

    done
    
fi

rm -f inputdir.txt
rm -f *.pcm *.d *.so

exit