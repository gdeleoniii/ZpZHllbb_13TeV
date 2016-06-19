#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_4_15/src
cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

ch=(ele mu)

for ((i=0; i<${#ch[@]}; i++)); do

    samplePath=/data7/htong/skim_NCUGlobalTuples

    cd $pwd/${ch[$i]}

    echo "We are now in " $PWD
    echo "Processing Z+jets background..."

    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\)
    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\)
    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\)
    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\)

    mv *root Zjets

    rm -f inputdir.txt
    rm -f *.pcm *.d *.so

    echo "Done. Move to next directory..."
    cd ../

done

echo "All the jobs are finished."

exit