#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

ch=(ele mu)
rg=(central up down)
samplePath=/data7/htong/skim_NCUGlobalTuples

for ((i=0; i<${#ch[@]}; i++)); do

    for ((j=0; j<${#rg[@]}; j++)); do

	if [ `echo ${ch[$i]} | grep -c "mu"` -gt 0 ]; then

	    echo "Processing muon data set..."

	    root -q -b -l jetEnergyTree.C\(\"$samplePath/skim_${ch[$i]}_SingleMuon_Run2015D-05Oct2015-v1.root\"\,\"SingleMuon-Run2015D-v1\"\,\"${rg[$j]}\"\,\"${ch[$i]}\"\)
	    root -q -b -l jetEnergyTree.C\(\"$samplePath/skim_${ch[$i]}_SingleMuon_Run2015D-PromptReco-v4.root\"\,\"SingleMuon-Run2015D-v4\"\,\"${rg[$j]}\"\,\"${ch[$i]}\"\)

	    mv *root data

	elif [ `echo ${ch[$i]} | grep -c "ele"` -gt 0 ]; then

	    echo "Processing electron data set..."

            root -q -b -l jetEnergyTree.C\(\"$samplePath/skim_${ch[$i]}_SingleElectron_Run2015D-05Oct2015-v1.root\"\,\"SingleElectron-Run2015D-v1\"\,\"${rg[$j]}\"\,\"${ch[$i]}\"\)
	    root -q -b -l jetEnergyTree.C\(\"$samplePath/skim_${ch[$i]}_SingleElectron_Run2015D-PromptReco-v4.root\"\,\"SingleElectron-Run2015D-v4\"\,\"${rg[$j]}\"\,\"${ch[$i]}\"\)

	    mv *root data

	fi

	echo "Processing VV background..."

	root -q -b -l jetEnergyTree.C\(\"$samplePath/skim_${ch[$i]}_crab_WW_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WW_TuneCUETP8M1_13TeV\"\,\"${rg[$j]}\"\,\"${ch[$i]}\"\)
	root -q -b -l jetEnergyTree.C\(\"$samplePath/skim_${ch[$i]}_crab_WZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WZ_TuneCUETP8M1_13TeV\"\,\"${rg[$j]}\"\,\"${ch[$i]}\"\)
	root -q -b -l jetEnergyTree.C\(\"$samplePath/skim_${ch[$i]}_crab_ZZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"ZZ_TuneCUETP8M1_13TeV\"\,\"${rg[$j]}\"\,\"${ch[$i]}\"\)

	echo "Processing ZH background..."

	root -q -b -l jetEnergyTree.C\(\"$samplePath/skim_${ch[$i]}_crab_ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8.root\"\,\"ZH_HToBB_ZToLL_M125_13TeV\"\,\"${rg[$j]}\"\,\"${ch[$i]}\"\)

	echo "Processing ttbar background..."

	root -q -b -l jetEnergyTree.C\(\"$samplePath/skim_${ch[$i]}_crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root\"\,\"TT_TuneCUETP8M1_13TeV\"\,\"${rg[$j]}\"\,\"${ch[$i]}\"\)

	mv *root minor

    done

    rm -f inputdir.txt
    rm -f *.pcm *.d *.so

    echo "Done."

done

echo "All the jobs are finished."

exit