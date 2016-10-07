#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src

## do cmsenv ##

cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

CHAN=(ele mu)
BTAG=(1 2)

## generate the necessary text file (contain event numbers) and root file (contain hist of ZH mass) ##

mass=(800 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000)

for ((i=0; i<${#CHAN[@]}; i++)); do
    for ((j=0; j<${#BTAG[@]}; j++)); do

	echo -e "*** Generate the necessary text file and root file ***"
	
	rootfile=mZH${CHAN[$i]}btag${BTAG[$j]}.root
	textfile=nEv${CHAN[$i]}btag${BTAG[$j]}.txt

	root -q -b -l nZHplots.C\(\"${CHAN[$i]}\"\,\"${BTAG[$j]}\"\,\"$rootfile\"\,\"$textfile\"\)
	
        ## make data cards for the combine tool ##
	
	dataCarddr=dataCards${CHAN[$i]}btag${BTAG[$j]}
	
	mkdir $dataCarddr
	
	echo -e "*** Make data cards for the combine tool by using: " $textfile " ***"
	echo -e "*** Data cards move to: " $dataCarddr " ***"
	
	uncertaintyfile=systUncOnSigEff/${CHAN[$i]}_${BTAG[$j]}btag_systUncOnSigEff.txt

	python MakeDataCards.py $textfile $rootfile ./$dataCarddr $uncertaintyfile
	
	rm -f DataCard_MXXXGeV.txt
	mv $rootfile $dataCarddr

        ## use the combine tool ##
	
	cd $cmsswdr/HiggsAnalysis/CombinedLimit/src
	
	for ((k=0; k<${#mass[@]}; k++)); do
	    
            dataCard=DataCard_M${mass[$k]}GeV_MonoHbb_13TeV.txt
            rootFile=higgsCombineCounting.Asymptotic.mZH${mass[$k]}.root
	    
            echo -e "*** Using data card: " $dataCard " to calculate limits ***"
	    
            combine -M Asymptotic $pwd/$dataCarddr/$dataCard
            mv higgsCombineTest.Asymptotic.mH120.root $pwd/$rootFile
	    
	done
	
        ## plot the results from the root files generate from combine tool ##
	
	cd $pwd
	echo -e "*** Plot the results using plotAsymptotic.C ***"
	
	root -q -b -l plotAsymptotic.C\(\"${CHAN[$i]}\"\,\"${BTAG[$j]}\"\)

	echo -e ""
	echo -e ""
	
    done
done

## combine data cards ##

combineCarddr=combineCards
mkdir $combineCarddr

eachCarddr=($(ls -d dataCards*))

for ((i=0; i<${#mass[@]}; i++)); do

    cd $pwd

    dataCard=DataCard_M${mass[$i]}GeV_MonoHbb_13TeV.txt
    combineCard=combine_DataCard_M${mass[$i]}GeV_MonoHbb_13TeV.txt

    echo -e "*** Combine data cards in " ${eachCarddr[0]} " " ${eachCarddr[1]} " " ${eachCarddr[2]} " " ${eachCarddr[3]} " ***"

    combineCards.py Name1=$pwd/${eachCarddr[0]}/$dataCard Name2=$pwd/${eachCarddr[1]}/$dataCard Name3=$pwd/${eachCarddr[2]}/$dataCard Name4=$pwd/${eachCarddr[3]}/$dataCard > $combineCard

    echo -e "*** Output card: " $combineCard " move to " $combineCarddr " ***"

    mv $combineCard $combineCarddr

    ## use the combine tool ##

    cd $cmsswdr/HiggsAnalysis/CombinedLimit/src

    rootFile=higgsCombineCounting.Asymptotic.combine.mZH${mass[$i]}.root

    echo -e "*** Using data card: " $combineCard " to calculate limits ***"

    combine -M Asymptotic $pwd/$combineCarddr/$combineCard
    mv higgsCombineTest.Asymptotic.mH120.root $pwd/$rootFile

done

## plot the results from the root files generate from combine tool ##

cd $pwd
echo -e "*** Plot the combine results using plotAsymptotic.C ***"

root -q -b -l plotAsymptotic.C++\(\"ele+mu\"\,\"1+2\"\)

## all jobs are completed ##

resultsdr=/afs/cern.ch/user/h/htong/www/limitResults

rm -rf $resultsdr
mkdir $resultsdr
mv *pdf $resultsdr
rm -f *.d *.so *.pcm 
rm -f higgsCombineCounting*root

echo -e "*** All jobs are completed ***"
echo -e ""

exit
