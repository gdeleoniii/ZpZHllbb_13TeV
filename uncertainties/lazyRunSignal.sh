#!/bin/sh

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.04.10/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

root -q -b -l btagging.signal/bTaggingUnc.C\(\)
root -q -b -l jetEnergyScale.signal/jetEnergyScale.C\(\)
root -q -b -l leptonScale.signal/leptonScaleUnc.C\(\)
root -q -b -l leptonTrigger.signal/leptonTriggerUnc.C\(\)
root -q -b -l pdfScale.signal/pdfScaleUnc.C\(\)
root -q -b -l pileup.signal/pileUpWeight.C\(\)

exit