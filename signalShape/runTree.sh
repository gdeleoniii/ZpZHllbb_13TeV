#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

ch=(ele mu)
mass=(800 1000 1200 1400 1600 1800 2000 2500 3000 3500)
samplePath=/data7/htong/skim_NCUGlobalTuples

for ((i=0; i<${#ch[@]}; i++)); do

    echo "Processing signal samples..."

    for mzh in "${mass[@]}"
    do

	signalName=skim_${ch[$i]}_crab_ZprimeToZhToZlephbb_narrow_M-${mzh}_13TeV-madgraph.root
	root -q -b -l signalTree.C\(\"$samplePath/$signalName\"\,\"ZprimeToZhToZlephbb_M-${mzh}_13TeV\"\,\"${ch[$i]}\"\,${mzh}\)

    done

    mv *root signal

    rm -f inputdir.txt
    rm -f *.pcm *.d *.so

    echo "Done."

done

echo "All the jobs are finished."

exit