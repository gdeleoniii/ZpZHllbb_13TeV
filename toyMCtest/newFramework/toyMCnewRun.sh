#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src
cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

channel=(ele mu)

for ((i=0; i<${#channel[@]}; i++)); do

    mcpath=/data7/htong/skim_samples/${channel[$i]}

    cd $pwd/${channel[$i]}

    echo "We are now in " $PWD
    echo "Processing DY+jets background..."

    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-100to200_13TeV\"\)
    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-200to400_13TeV\"\)
    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-400to600_13TeV\"\)
    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\"\,\"DYJetsToLL_M-50_HT-600toInf_13TeV\"\)

    mv *root DYjets

    echo "Processing diBosons background..."

    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_WW_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WW_TuneCUETP8M1_13TeV\"\)
    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_WZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"WZ_TuneCUETP8M1_13TeV\"\)
    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_ZZ_TuneCUETP8M1_13TeV-pythia8.root\"\,\"ZZ_TuneCUETP8M1_13TeV\"\)

    mv *root diBosons

    echo "Processing ttbar background..."

    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root\"\,\"TT_TuneCUETP8M1_13TeV\"\)

    mv *root ttbar

    echo "Processing singleTop background..."

    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root\"\,\"ST_s_4f_leptonDecays_13TeV\"\)
    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root\"\,\"ST_t_antitop_4f_leptonDecays_13TeV\"\)
    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root\"\,\"ST_t_top_4f_leptonDecays_13TeV\"\)
    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root\"\,\"ST_tW_antitop_5f_inclusiveDecays_13TeV\"\)
    root -q -b -l toyMCnew_${channel[$i]}.C+\(\"$mcpath/skim_${channel[$i]}_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root\"\,\"ST_tW_top_5f_inclusiveDecays_13TeV\"\)

    mv *root singleTop

    rm -f inputdir.txt
    rm -f *.pcm *.d *.so

    echo "Done. Move to next directory..."
    cd ../

done

echo "All the jobs are finished."

exit