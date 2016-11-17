#!/bin/sh

if [ -z $1 ]; then
    echo "Usage: $0 [muon or electron]"
    exit 0
fi

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

ana=/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/skimSamples
dataPath=/data7/syu/NCUGlobalTuples
mcPath=/data7/htong/NCUGlobalTuples
chan=$1

tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$dataPath/Run2015D/b6fa618/SingleElectron_Run2015D-05Oct2015-v1/\"\,0\)  &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$dataPath/Run2015D/b6fa618/SingleElectron_Run2015D-PromptReco-v4/\"\,0\) &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$dataPath/Run2015D/b6fa618/SingleMuon_Run2015D-05Oct2015-v1/\"\,0\)      &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$dataPath/Run2015D/b6fa618/SingleMuon_Run2015D-PromptReco-v4/\"\,0\)     &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/\"\,0\) &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/\"\,0\) &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/\"\,0\) &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/\"\,0\) &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8/\"\,0\)                &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_WW_TuneCUETP8M1_13TeV-pythia8/\"\,0\)                       &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_WZ_TuneCUETP8M1_13TeV-pythia8/\"\,0\)                       &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZZ_TuneCUETP8M1_13TeV-pythia8/\"\,0\)                       &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/\"\,0\)            &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-600_13TeV-madgraph/\"\,0\)     &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-800_13TeV-madgraph/\"\,0\)     &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-1000_13TeV-madgraph/\"\,0\)    &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-1200_13TeV-madgraph/\"\,0\)    &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-1400_13TeV-madgraph/\"\,0\)    &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-1600_13TeV-madgraph/\"\,0\)    &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-1800_13TeV-madgraph/\"\,0\)    &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-2000_13TeV-madgraph/\"\,0\)    &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-2500_13TeV-madgraph/\"\,0\)    &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-3000_13TeV-madgraph/\"\,0\)    &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-3500_13TeV-madgraph/\"\,0\)    &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-4000_13TeV-madgraph/\"\,0\)    &> Log &
tmp=`mktemp -d --tmpdir=.` && cd $tmp && cp $ana/* . && root -q -b -l runSkimTree.C\(\"$chan\"\,\"$mcPath/crab_ZprimeToZhToZlephbb_narrow_M-4500_13TeV-madgraph/\"\,0\)    &> Log &

wait

mv tmp.*/skim*root .

exit