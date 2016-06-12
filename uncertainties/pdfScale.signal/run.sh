#!/bin/sh

pwd=$PWD
cmsswdr=/afs/cern.ch/work/h/htong/CMSSW_7_4_15/src
cd $cmsswdr
export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`
cd $pwd

root -q -b -l pdfScaleUnc.C+\(1\,0\,2\,1\,\"mur1\"\)
root -q -b -l pdfScaleUnc.C+\(2\,0\,2\,1\,\"mur1\"\)
root -q -b -l pdfScaleUnc.C+\(1\,0\,6\,2\,\"muf1\"\)
root -q -b -l pdfScaleUnc.C+\(2\,0\,6\,2\,\"muf1\"\)
root -q -b -l pdfScaleUnc.C+\(1\,10\,109\,1\,\"pdf\"\)
root -q -b -l pdfScaleUnc.C+\(2\,10\,109\,1\,\"pdf\"\)

dr=pdfScaleResults

rm -f *.pcm *.d *.so
rm -rf $HOME/www/$dr
mkdir $dr
mv *txt $dr
mv $dr $HOME/www

echo "All the jobs are finished."

exit