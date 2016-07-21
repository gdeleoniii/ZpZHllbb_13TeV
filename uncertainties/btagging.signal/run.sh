#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

root -q -b -l bTagUnc.C\(\)

mkdir signalBtagResults/
mv *txt signalBtagResults/
rm -rf $HOME/www/signalBtagResults/
mv signalBtagResults/ $HOME/www/

exit
