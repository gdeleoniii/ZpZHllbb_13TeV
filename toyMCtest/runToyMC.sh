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

    if [ `echo ${ch[$i]} | grep -c "mu"` -gt 0 ]; then

	echo "Processing muon data set..."

	root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_SingleMuon_Run2015D-05Oct2015-v1.root\"\,\"SingleMuon-Run2015D-v1\"\)
	root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_SingleMuon_Run2015D-PromptReco-v4.root\"\,\"SingleMuon-Run2015D-v4\"\)

	mv *root data

    elif [ `echo ${ch[$i]} | grep -c "ele"` -gt 0 ]; then

	echo "Processing electron data set..."

        root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_SingleElectron_Run2015D-05Oct2015-v1.root\"\,\"SingleElectron-Run2015D-v1\"\)
	root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_SingleElectron_Run2015D-PromptReco-v4.root\"\,\"SingleElectron-Run2015D-v4\"\)

	mv *root data

    fi

    echo "Processing Z+jets background..."

    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\)
    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\)
    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\)
    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\)

    mv *root Zjets

    echo "Processing VV background..."

    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_WW_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WW_TuneCUETP8M1_13TeV\"\)
    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_WZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WZ_TuneCUETP8M1_13TeV\"\)
    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_ZZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"ZZ_TuneCUETP8M1_13TeV\"\)

    mv *root VV

    echo "Processing ZH background..."

    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8.root\"\,\"ZH_HToBB_ZToLL_M125_13TeV\"\)

    mv *root ZH

    echo "Processing ttbar background..."

    root -q -b -l toyMC_${ch[$i]}.C+\(\"$samplePath/skim_${ch[$i]}_crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root\"\,\"TT_TuneCUETP8M1_13TeV\"\)

    mv *root TT

    rm -f inputdir.txt
    rm -f *.pcm *.d *.so

    echo "Done. Move to next directory..."
    cd ../

done

echo "All the jobs are finished."

exit