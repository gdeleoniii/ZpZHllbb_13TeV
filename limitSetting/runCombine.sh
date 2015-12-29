#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src

## do cmsenv ##

cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

## generate the necessary text file (contain event numbers) and root file (contain hist of ZH mass) ##

mkdir outputmZHSig/catOne; cp outputmZHSig/*catOne.root outputmZHSig/catOne
mkdir outputmZHSig/catTwo; cp outputmZHSig/*catTwo.root outputmZHSig/catTwo

# catOne #

catOnetextfile=nEventZHcatOne.txt
catOnerootfile=mZHmuSetLimitcatOne.root

echo -e "*** Generate the necessary text file and root file ***"
root -q -b -l mZHSIGplots.C++\(\"outputmZHSig/catOne\"\,\"$catOnerootfile\"\,\"$catOnetextfile\"\)

# catTwo #

catTwotextfile=nEventZHcatTwo.txt
catTworootfile=mZHmuSetLimitcatTwo.root

echo -e "*** Generate the necessary text file and root file ***"
root -q -b -l mZHSIGplots.C++\(\"outputmZHSig/catTwo\"\,\"$catTworootfile\"\,\"$catTwotextfile\"\)

## check are the necessary files exist ##

if [ -e $catOnetextfile ] && [ -e $catTwotextfile ] && [ -e $catOnerootfile ] && [ -e $catTworootfile ]
then
    echo -e "*** The necessary text file " $catOnetextfile " and " $catTwotextfile " are exist! ***"
    echo -e "*** The necessary root file " $catOnerootfile " and " $catTworootfile " are exist! ***"
else
    echo -e "*** The necessary text file " $catOnetextfile " and " $catTwotextfile " doesn't exist! ***"
    echo -e "*** The necessary root file " $catOnerootfile " and " $catTworootfile " doesn't exist! ***"
    exit 0
fi

## make data cards for the combine tool ##

# catOne cards #

catOnecarddr=catOnedataCards

if [ -d $catOnecarddr ]
then
    echo -e "*** Data card directory is " $catOnecarddr " ***"
else
    echo -e "*** Generate data card directory: " $catOnecarddr " ***"
    mkdir $catOnecarddr
fi

echo -e "*** Make data cards for the combine tool by using: " $catOnetextfile " ***"
echo -e "*** Data cards move to: " $catOnecarddr " ***"

python MakeDataCards.py $catOnetextfile $catOnerootfile ./$catOnecarddr
rm -f DataCard_MXXXGeV.txt
mv $catOnerootfile $catOnecarddr

# catTwo cards #

catTwocarddr=catTwodataCards

if [ -d $catTwocarddr ]
then
    echo -e "Data card directory is " $catTwocarddr
else
    echo -e "Generate data card directory: " $catTwocarddr
    mkdir $catTwocarddr
fi

echo -e "*** Make data cards for the combine tool by using: " $catTwotextfile " ***"
echo -e "*** Data cards move to: " $catTwocarddr " ***"

python MakeDataCards.py $catTwotextfile $catTworootfile ./$catTwocarddr
rm -f DataCard_MXXXGeV.txt
mv $catTworootfile $catTwocarddr

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

    echo -e "*** Combine cards: " $catTwocarddr"/"$dataCard " and " $catOnecarddr"/"$dataCard " ***"

    combineCards.py Name1=$pwd/$catTwocarddr/$dataCard Name2=$pwd/$catOnecarddr/$dataCard > $combineCard

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
rm -rf outputmZHSig/catOne outputmZHSig/catTwo
rm -f higgsCombineCounting*root

echo -e "*** All jobs are completed ***"

exit
