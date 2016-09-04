#!/bin/bash

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

root -q -b -l btagging.signal/bTaggingUnc.C\(\)
root -q -b -l jetEnergyScale.signal/jetEnergyScale.C\(\)
root -q -b -l leptonScale.signal/leptonScaleUnc.C\(\)
root -q -b -l leptonTrigger.signal/leptonTriggerUnc.C\(\)
root -q -b -l pdfScale.signal/pdfScaleUnc.C\(\)
root -q -b -l pileup.signal/pileUpWeight.C\(\)

# Below part is to combine all the results into one text file

targetDir=signal*Results
resultsDir=systUncOnSigEff
channel=(ele mu)
cat=(1 2)

for ((i=0; i<${#channel[@]}; i++)); do
    for ((j=0; j<${#cat[@]}; j++)); do

	paste $targetDir/${channel[$i]}_${cat[$j]}*Unc.txt | awk '{{printf("%s\t%s\t",$1,$2)};k=3;while($k){printf("%s\t",$k);k+=3}printf("\n")}' > ${channel[$i]}_${cat[$j]}btag_systUncOnSigEff.txt

    done
done

rm -rf $targetDir
mkdir -p $resultsDir
mv *systUncOnSigEff.txt $resultsDir

exit