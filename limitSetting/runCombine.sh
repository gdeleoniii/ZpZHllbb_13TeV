#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src

## do cmsenv ##

cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

## generate the necessary text file (contain event numbers) and root file (contain hist of ZH mass) ##

# catEle #

catEletextfile=nEventZHcatEle.txt
catElerootfile=mZHmuSetLimitcatEle.root

echo -e "*** Generate the necessary text file and root file ***"
root -q -b -l nZHplots.C++\(\"Electron\"\,\"ele\"\,\"$catElerootfile\"\,\"$catEletextfile\"\)

# catMu #

catMutextfile=nEventZHcatMu.txt
catMurootfile=mZHmuSetLimitcatMu.root

echo -e "*** Generate the necessary text file and root file ***"
root -q -b -l nZHplots.C++\(\"Muon\"\,\"mu\"\,\"$catMurootfile\"\,\"$catMutextfile\"\)

## check are the necessary files exist ##

if [ -e $catEletextfile ] && [ -e $catMutextfile ] && [ -e $catElerootfile ] && [ -e $catMurootfile ]
then
    echo -e "*** The necessary text file " $catEletextfile " and " $catMutextfile " are exist! ***"
    echo -e "*** The necessary root file " $catElerootfile " and " $catMurootfile " are exist! ***"
else
    echo -e "*** The necessary text file " $catEletextfile " and " $catMutextfile " doesn't exist! ***"
    echo -e "*** The necessary root file " $catElerootfile " and " $catMurootfile " doesn't exist! ***"
    exit 0
fi

## make data cards for the combine tool ##

# catEle cards #

catEleCarddr=catEledataCards

if [ -d $catEleCarddr ]
then
    echo -e "*** Data card directory is " $catEleCarddr " ***"
else
    echo -e "*** Generate data card directory: " $catEleCarddr " ***"
    mkdir $catEleCarddr
fi

echo -e "*** Make data cards for the combine tool by using: " $catEletextfile " ***"
echo -e "*** Data cards move to: " $catEleCarddr " ***"

python MakeDataCards.py $catEletextfile $catElerootfile ./$catEleCarddr
rm -f DataCard_MXXXGeV.txt
mv $catElerootfile $catEleCarddr

# catMu cards #

catMuCarddr=catMudataCards

if [ -d $catMuCarddr ]
then
    echo -e "Data card directory is " $catMuCarddr
else
    echo -e "Generate data card directory: " $catMuCarddr
    mkdir $catMuCarddr
fi

echo -e "*** Make data cards for the combine tool by using: " $catMutextfile " ***"
echo -e "*** Data cards move to: " $catMuCarddr " ***"

python MakeDataCards.py $catMutextfile $catMurootfile ./$catMuCarddr
rm -f DataCard_MXXXGeV.txt
mv $catMurootfile $catMuCarddr

## combine data cards ##

dataCarddr=dataCards

if [ -d $dataCarddr ]
then
    echo -e "Combine data card directory is " $dataCarddr
else
    echo -e "Generate data card directory: " $dataCarddr
    mkdir $dataCarddr
fi

massPoints=(800 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000)

for ((i=0; i<${#massPoints[@]}; i++)); do

    dataCard=DataCard_M${massPoints[$i]}GeV_MonoHbb_13TeV.txt
    combineCard=combine_DataCard_M${massPoints[$i]}GeV_MonoHbb_13TeV.txt

    echo -e "*** Combine cards: " $catMuCarddr"/"$dataCard " and " $catEleCarddr"/"$dataCard " ***"

    combineCards.py Name1=$pwd/$catMuCarddr/$dataCard Name2=$pwd/$catEleCarddr/$dataCard > $combineCard

    echo -e "*** Output card: " $combineCard " move to " $dataCarddr " ***"

    mv $combineCard $dataCarddr

done

## use the combine tool ##

cd $cmsswdr/HiggsAnalysis/CombinedLimit/src

for ((i=0; i<${#massPoints[@]}; i++)); do

    dataCard=combine_DataCard_M${massPoints[$i]}GeV_MonoHbb_13TeV.txt
    rootFile=higgsCombineCounting.Asymptotic.mZH${massPoints[$i]}.root

    echo -e "*** Using data card: " $dataCard " to calculate limits ***"

    combine -M Asymptotic $pwd/$dataCarddr/$dataCard
    mv higgsCombineTest.Asymptotic.mH120.root $pwd/$rootFile

done

## plot the results from the root files generate from combine tool ##

cd $pwd
echo -e "*** Plot the results using plotAsymptotic.C ***"

root -q -b -l plotAsymptotic.C++\(\)

## all jobs are completed ##

mv *pdf /afs/cern.ch/user/h/htong/www
rm -f *.d *.so *.pcm 
rm -rf $catEleCarddr $catMuCarddr
rm -f higgsCombineCounting*root

echo -e "*** All jobs are completed ***"

exit
