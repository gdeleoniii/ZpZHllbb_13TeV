#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

pwd=$PWD

ch=(ele mu)
rg=(central up down)

for ((i=0; i<${#ch[@]}; i++)); do

    samplePath=/data7/htong/skim_NCUGlobalTuples

    cd $pwd/${ch[$i]}

    echo "We are now in " $PWD
    echo "Processing Z+jets background..."

    for ((j=0; j<${#rg[@]}; j++)); do

	root -q -b -l toyMC_${ch[$i]}.C\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\,\"${rg[$j]}\"\)
	root -q -b -l toyMC_${ch[$i]}.C\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\,\"${rg[$j]}\"\)
	root -q -b -l toyMC_${ch[$i]}.C\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\,\"${rg[$j]}\"\)
	root -q -b -l toyMC_${ch[$i]}.C\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\,\"${rg[$j]}\"\)

    done
    
    mv *root Zjets
    
    rm -f inputdir.txt
    rm -f *.pcm *.d *.so
    
    echo "Done. Move to next directory..."
    cd ../
    
done

echo "All the jobs are finished."

exit