#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src
cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

ch=(ele mu)

for ((i=0; i<${#ch[@]}; i++)); do

    mcpath=/data7/htong/skim_samples/${ch[$i]}
    datapath=/data7/htong/

    cd $pwd/${ch[$i]}

    echo "We are now in " $PWD

    if [ `echo ${ch[$i]} | grep -c "mu"` -gt 0 ]; then

	echo "Processing muon data set..."

	root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_crab_SingleMuon-Run2015D-05Oct2015-v1_20151119_2p2fb_SingleMuTextFile.root\"\,\"SingleMuon-Run2015D-05Oct2015-v1\"\)
	root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_crab_SingleMuon-Run2015D-PromptReco-V420151119_2p2fb_SingleMuTextFile.root\"\,\"SingleMuon-Run2015D-PromptReco-V4\"\)

	mv *root data

    elif [ `echo ${ch[$i]} | grep -c "ele"` -gt 0 ]; then

	echo "Processing electron data set..."

        root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_crab_SingleElectron-Run2015D-05Oct2015-v1_20151117_2p2fb_SingleEleTextFile.root\"\,\"SingleElectron-Run2015D-05Oct2015-v1\"\)
	root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_crab_SingleElectron-Run2015D-PromptReco-V420151117_2p2fb_SingleEleTextFile.root\"\,\"SingleElectron-Run2015D-PromptReco-V4\"\)

	mv *root data

    fi

    echo "Processing DY+jets background..."

    root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\)
    root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\)
    root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\)
    root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\)

    mv *root Zjets

    echo "Processing diBosons background..."

    root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_WW_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WW_TuneCUETP8M1_13TeV\"\)
    root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_WZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WZ_TuneCUETP8M1_13TeV\"\)
    root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_ZZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"ZZ_TuneCUETP8M1_13TeV\"\)

    mv *root VV

    echo "Processing ttbar background..."

    root -q -b -l toyMCnew_${ch[$i]}.C+\(\"$mcpath/skim_${ch[$i]}_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root\"\,\"TT_TuneCUETP8M1_13TeV\"\)

    mv *root TT

    rm -f inputdir.txt
    rm -f *.pcm *.d *.so

    echo "Done. Move to next directory..."
    cd ../

done

echo "All the jobs are finished."

exit