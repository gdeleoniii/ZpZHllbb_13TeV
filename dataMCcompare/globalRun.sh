#!/bin/sh

if [ -z $2 ]; then

    echo "Usage: $0 [macro without .C] [channel]"
    exit 0

fi

channel=$2

#Do cmsenv:
#pwd=$PWD
#cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src
#cd $cmsswdr
#export SCRAM_ARCH=slc6_amd64_gcc481
#eval `scramv1 runtime -sh`
#cd $pwd

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

samplePath=/data7/htong/skim_NCUGlobalTuples

if [ $channel == ele ]; then

    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_SingleElectron_Run2015D-05Oct2015-v1.root\"\,\"SingleElectron-Run2015D-v1\"\)
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_SingleElectron_Run2015D-PromptReco-v4.root\"\,\"SingleElectron-Run2015D-v4\"\)

elif [ $channel == mu ]; then

    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_SingleMuon_Run2015D-05Oct2015-v1.root\"\,\"SingleMuon-Run2015D-v1\"\)             
    root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_SingleMuon_Run2015D-PromptReco-v4.root\"\,\"SingleMuon-Run2015D-v4\"\)

fi

root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\)
root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\)
root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\)
root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\)
root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_WW_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WW_TuneCUETP8M1_13TeV\"\)
root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_WZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WZ_TuneCUETP8M1_13TeV\"\)
root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_ZZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"ZZ_TuneCUETP8M1_13TeV\"\)
root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8.root\"\,\"ZH_HToBB_ZToLL_M125_13TeV\"\)
root -q -b -l $1.C+\(\"$samplePath/skim_${channel}_crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root\"\,\"TT_TuneCUETP8M1_13TeV\"\)

rm -f inputdir.txt
rm -f *.pcm *.d *.so

exit