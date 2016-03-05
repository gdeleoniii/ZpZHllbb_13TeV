#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_1_5/src
cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

channel=(mu ele)
cat=(1 2)
phy=(DYjets TTbar Dibosons SingleTop)

for ((i=0; i<${#channel[@]}; i++)); do
    for ((j=0; j<${#cat[@]}; j++)); do
	for ((k=0; k<${#phy[@]}; k++)); do
	    root -q -b -l for${phy[$k]}.C+\(\"${channel[$i]}/cat${cat[$j]}\"\,\"${channel[$i]}_cat${cat[$j]}_${phy[$k]}\"\)
	done
    done
done

mv *pdf fitTestResults/
echo "All the jobs are finished."

exit