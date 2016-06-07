#!/bin/sh

if [ -z $2 ]; then

    echo "Usage: $0 [channel] [btag]"
    exit 0

fi

chan=$1
btag=$2

samplePath=/data7/htong/skim_NCUGlobalTuples

if [ $chan == ele ]; then

    root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_SingleElectron_Run2015D-05Oct2015-v1.root\"\,\"SingleElectron-Run2015D-v1\"\,$btag\)
    root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_SingleElectron_Run2015D-PromptReco-v4.root\"\,\"SingleElectron-Run2015D-v4\"\,$btag\)

elif [ $chan == mu ]; then

    root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_SingleMuon_Run2015D-05Oct2015-v1.root\"\,\"SingleMuon-Run2015D-v1\"\,$btag\)             
    root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_SingleMuon_Run2015D-PromptReco-v4.root\"\,\"SingleMuon-Run2015D-v4\"\,$btag\)

fi

mass=(800 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000)

for i in "${mass[@]}"
do
    signalName=skim_${chan}_crab_ZprimeToZhToZlephbb_narrow_M-${i}_13TeV-madgraph.root
    root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/$signalName\"\,\"ZprimeToZhToZlephbb_M-${i}_13TeV\"\,$btag\)
done

root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\,$btag\)
root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_crab_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\,$btag\)
root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\,$btag\)
root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_crab_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\,$btag\)
root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_crab_WW_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WW_TuneCUETP8M1_13TeV\"\,$btag\)
root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_crab_WZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WZ_TuneCUETP8M1_13TeV\"\,$btag\)
root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_crab_ZZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"ZZ_TuneCUETP8M1_13TeV\"\,$btag\)
root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_crab_ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8.root\"\,\"ZH_HToBB_ZToLL_M125_13TeV\"\,$btag\)
root -q -b -l mZH${chan}Limit.C+\(\"$samplePath/skim_${chan}_crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root\"\,\"TT_TuneCUETP8M1_13TeV\"\,$btag\)

outputName=output_${chan}_${btag}btag
mkdir $outputName
mv *_mZHLimit.root $outputName

rm -f inputdir.txt
rm -f *.pcm *.d *.so

exit