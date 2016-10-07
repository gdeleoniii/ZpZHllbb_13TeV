#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

chan=(ele mu)
btag=(1 2)
mass=(800 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000)
samplePath=/data7/htong/skim_NCUGlobalTuples

for ((i=0; i<${#chan[@]}; i++)); do
    for ((j=0; j<${#btag[@]}; j++)); do

	echo "Now running ==> channel: " ${chan[$i]} " category: " ${btag[$j]}

	if [ ${chan[$i]} == ele ]; then

	    root -q -b -l mZHLimit.C\(\"$samplePath/skim_ele_SingleElectron_Run2015D-05Oct2015-v1.root\"\,\"SingleElectron-Run2015D-v1\"\,\"ele\"\,${btag[$j]}\)
	    root -q -b -l mZHLimit.C\(\"$samplePath/skim_ele_SingleElectron_Run2015D-PromptReco-v4.root\"\,\"SingleElectron-Run2015D-v4\"\,\"ele\"\,${btag[$j]}\)

	elif [ ${chan[$i]} == mu ]; then

	    root -q -b -l mZHLimit.C\(\"$samplePath/skim_mu_SingleMuon_Run2015D-05Oct2015-v1.root\"\,\"SingleMuon-Run2015D-v1\"\,\"mu\"\,${btag[$j]}\)
	    root -q -b -l mZHLimit.C\(\"$samplePath/skim_mu_SingleMuon_Run2015D-PromptReco-v4.root\"\,\"SingleMuon-Run2015D-v4\"\,\"mu\"\,${btag[$j]}\)

	fi
	
	for mzh in "${mass[@]}"; do
	    signalName=skim_${chan[$i]}_crab_ZprimeToZhToZlephbb_narrow_M-${mzh}_13TeV-madgraph.root
	    root -q -b -l mZHLimit.C\(\"$samplePath/$signalName\"\,\"ZprimeToZhToZlephbb_M-${mzh}_13TeV\"\,\"${chan[$i]}\"\,${btag[$j]}\,${mzh}\)
	done

	root -q -b -l mZHLimit.C\(\"$samplePath/skim_${chan[$i]}_crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\,\"${chan[$i]}\"\,${btag[$j]}\)
	root -q -b -l mZHLimit.C\(\"$samplePath/skim_${chan[$i]}_crab_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\,\"${chan[$i]}\"\,${btag[$j]}\)
	root -q -b -l mZHLimit.C\(\"$samplePath/skim_${chan[$i]}_crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\,\"${chan[$i]}\"\,${btag[$j]}\)
	root -q -b -l mZHLimit.C\(\"$samplePath/skim_${chan[$i]}_crab_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\,\"${chan[$i]}\"\,${btag[$j]}\)
	root -q -b -l mZHLimit.C\(\"$samplePath/skim_${chan[$i]}_crab_WW_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WW_TuneCUETP8M1_13TeV\"\,\"${chan[$i]}\"\,${btag[$j]}\)
	root -q -b -l mZHLimit.C\(\"$samplePath/skim_${chan[$i]}_crab_WZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WZ_TuneCUETP8M1_13TeV\"\,\"${chan[$i]}\"\,${btag[$j]}\)
	root -q -b -l mZHLimit.C\(\"$samplePath/skim_${chan[$i]}_crab_ZZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"ZZ_TuneCUETP8M1_13TeV\"\,\"${chan[$i]}\"\,${btag[$j]}\)
	root -q -b -l mZHLimit.C\(\"$samplePath/skim_${chan[$i]}_crab_ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8.root\"\,\"ZH_HToBB_ZToLL_M125_13TeV\"\,\"${chan[$i]}\"\,${btag[$j]}\)
	root -q -b -l mZHLimit.C\(\"$samplePath/skim_${chan[$i]}_crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root\"\,\"TT_TuneCUETP8M1_13TeV\"\,\"${chan[$i]}\"\,${btag[$j]}\)

	outputName=output_${chan[$i]}_${btag[$j]}btag
	mkdir $outputName
	mv *_mZHLimit.root $outputName
	
    done
done

rm -f inputdir.txt
rm -f *.pcm *.d *.so

echo "All the jobs are finished."

exit