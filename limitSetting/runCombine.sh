#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src

## do cmsenv ##

cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

## generate the necessary text file (contain event numbers) and root file (contain hist of ZH mass) ##

mkdir outputmZHSig/loose; cp outputmZHSig/*loose.root outputmZHSig/loose
mkdir outputmZHSig/tight; cp outputmZHSig/*tight.root outputmZHSig/tight

# loose #

loosetextfile=nEventZHloose.txt
looserootfile=mZHmuSetLimitloose.root

echo -e "*** Generate the necessary text file and root file ***"
root -q -b -l mZHSIGplots.C++\(\"outputmZHSig/loose\"\,\"$looserootfile\"\,\"$loosetextfile\"\)

# tight #

tighttextfile=nEventZHtight.txt
tightrootfile=mZHmuSetLimittight.root

echo -e "*** Generate the necessary text file and root file ***"
root -q -b -l mZHSIGplots.C++\(\"outputmZHSig/tight\"\,\"$tightrootfile\"\,\"$tighttextfile\"\)

## check are the necessary files exist ##

if [ -e $loosetextfile ] && [ -e $tighttextfile ] && [ -e $looserootfile ] && [ -e $tightrootfile ]
then
    echo -e "*** The necessary text file " $loosetextfile " and " $tighttextfile " are exist! ***"
    echo -e "*** The necessary root file " $looserootfile " and " $tightrootfile " are exist! ***"
else
    echo -e "*** The necessary text file " $loosetextfile " and " $tighttextfile " doesn't exist! ***"
    echo -e "*** The necessary root file " $looserootfile " and " $tightrootfile " doesn't exist! ***"
    exit 0
fi

## make data cards for the combine tool ##

# loose cards #

loosecarddr=loosedataCards

if [ -d $loosecarddr ]
then
    echo -e "*** Data card directory is " $loosecarddr " ***"
else
    echo -e "*** Generate data card directory: " $loosecarddr " ***"
    mkdir $loosecarddr
fi

echo -e "*** Make data cards for the combine tool by using: " $loosetextfile " ***"
echo -e "*** Data cards move to: " $loosecarddr " ***"

python MakeDataCards.py $loosetextfile $looserootfile ./$loosecarddr
rm -f DataCard_MXXXGeV.txt
mv $looserootfile $loosecarddr

# tight cards #

tightcarddr=tightdataCards

if [ -d $tightcarddr ]
then
    echo -e "Data card directory is " $tightcarddr
else
    echo -e "Generate data card directory: " $tightcarddr
    mkdir $tightcarddr
fi

echo -e "*** Make data cards for the combine tool by using: " $tighttextfile " ***"
echo -e "*** Data cards move to: " $tightcarddr " ***"

python MakeDataCards.py $tighttextfile $tightrootfile ./$tightcarddr
rm -f DataCard_MXXXGeV.txt
mv $tightrootfile $tightcarddr

## combine data cards ##

datacarddr=dataCards

if [ -d $datacarddr ]
then
    echo -e "Combine data card directory is " $datacarddr
else
    echo -e "Generate data card directory: " $datacarddr
    mkdir $datacarddr
fi

massPoints=(800 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000)

for ((i=0; i<${#massPoints[@]}; i++)); do

    dataCard=DataCard_M${massPoints[$i]}GeV_MonoHbb_13TeV.txt
    combineCard=combine_DataCard_M${massPoints[$i]}GeV_MonoHbb_13TeV.txt

    echo -e "*** Combine cards: " $tightcarddr"/"$dataCard " and " $loosecarddr"/"$dataCard " ***"

    combineCards.py Name1=$pwd/$tightcarddr/$dataCard Name2=$pwd/$loosecarddr/$dataCard > $combineCard

    echo -e "*** Output card: " $combineCard " move to " $datacarddr " ***"

    mv $combineCard $datacarddr

done

## use the combine tool ##

cd $cmsswdr/HiggsAnalysis/CombinedLimit/src

for ((i=0; i<${#massPoints[@]}; i++)); do

    dataCard=combine_DataCard_M${massPoints[$i]}GeV_MonoHbb_13TeV.txt
    rootFile=higgsCombineCounting.Asymptotic.mZH${massPoints[$i]}.root

    echo -e "*** Using data card: " $dataCard " to calculate limits ***"

    combine -M Asymptotic $pwd/$datacarddr/$dataCard
    mv higgsCombineTest.Asymptotic.mH120.root $pwd/$rootFile

done

## plot the results from the root files generate from combine tool ##

cd $pwd
echo -e "*** Plot the results using plot_Asymptotic.C ***"

root -q -b -l plot_Asymptotic.C++\(\"Counting\"\)

## all jobs are completed ##

mv *pdf /afs/cern.ch/user/h/htong/www
rm -f *.d *.so *.pcm
rm -rf outputmZHSig/loose outputmZHSig/tight
rm -f higgsCombineCounting*root

echo -e "*** All jobs are completed ***"

exit
