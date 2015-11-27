#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src

## do cmsenv

cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

## generate the necessary text file (contain event numbers) and root file (contain hist of ZH mass)

textfile=nEventZH.txt
rootfile=mZHmuSetLimit.root

root -q -b -l mZHSIGplots.C++\(\"../muon/forLimitCalc/outputmZHmu/signalloose\"\,\"$rootfile\"\,\"$textfile\"\)

## check is the files exist

datacarddr=dataCards

if [ -d $datacarddr ]
then
    echo "Data card directory is " $datacarddr
else
    echo "Generate data card directory: " $datacarddr
    mkdir $datacarddr
fi

if [ -e $textfile ] && [ -e $rootfile ]
then
    echo "The necessary text file " $textfile " and root file " $rootfile " are exist!"
else
    echo "The necessary text file " $textfile " and root file " $rootfile " doesn't exist!"
    exit 0
fi

## make data cards for the combine tool

python MakeDataCards.py $textfile $rootfile ./$datacarddr
rm -f DataCard_MXXXGeV.txt
mv $rootfile $datacarddr

## use the combine tool

cd $cmsswdr/HiggsAnalysis/CombinedLimit/src

massPoints=(800 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000)

for ((i=0; i<${#massPoints[@]}; i++)); do

    dataCard=DataCard_M${massPoints[$i]}GeV_MonoHbb_13TeV.txt
    rootFile=higgsCombineCounting.Asymptotic.mZH${massPoints[$i]}.root

    echo -e "\n****** Using data card: " $dataCard " ******\n"

    combine -M Asymptotic $pwd/$datacarddr/$dataCard
    mv higgsCombineTest.Asymptotic.mH120.root $pwd/$rootFile

done

## plot the results from the root files generate from combine tool 

cd $pwd
root -q -b -l plot_Asymptotic.C++\(\"Counting\"\)

## all jobs are completed

mv *pdf /afs/cern.ch/user/h/htong/www
rm -f *.d *.so *.pcm

exit
