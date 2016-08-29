import sys
import os
#import ROOT
#ROOT.gROOT.SetBatch(True)  
#from DataSetInfo import *

if len(sys.argv) < 4 :
    print "insufficient options provided see help function "
    exit (1)

if len(sys.argv) == 4 :
    print ('You are making datacards for '+sys.argv[1]+' with '+sys.argv[2]+' and datacards will be saved in '+sys.argv[3])

inputtextfilename=sys.argv[1]
inputrootfilename=sys.argv[2]
dirtosave=sys.argv[3]

os.system('mkdir -p '+dirtosave)

## prepare template datacard and store in this variable. 
## do whatever change you want to do 
## and keep variables in capital letters with prefix "T" e.g. signal strength can be written as TSIGNAL and DYJets can be TDYJETS
## these can be replaced by either a python or shell script to make datacards for a specific mass point and analysis. 
templatecard='''
imax    1        number of channels
jmax    *        number of backgrounds
kmax    *        number of nuisance parameters (sources of systematical uncertainties)

-----------------------------------------------------------------------------------------------------------------

shapes  *        ZPTOZHMM  INPUTROOTFILE  $PROCESS  $PROCESS_$SYSTEMATIC

-----------------------------------------------------------------------------------------------------------------

bin              ZPTOZHMM
observation      DATARATE

-----------------------------------------------------------------------------------------------------------------

bin                          ZPTOZHMM     ZPTOZHMM     ZPTOZHMM    ZPTOZHMM    ZPTOZHMM    ZPTOZHMM   ZPTOZHMM
process                      SIGNAL       DYJETS       TTBAR       WW          WZ          ZZ         ZH
process                      0            1            2           3           4           5          6
rate                         SIGNALRATE   DYJETSRATE   TTBARRATE   WWRATE      WZRATE      ZZRATE     ZHRATE

-----------------------------------------------------------------------------------------------------------------

Alpha_bTag      shape        -          1.07       -          -          -          -          -   
Alpha_QCD       shape        -          1.07       -          -          -          -          -   
Alpha_PDF       shape        -          1.07       -          -          -          -          -   
Alpha_JES       shape        -          1.07       -          -          -          -          -   
SigEf_bTag      lnN          1.10       -          -          -          -          -          -
SigEf_QCD       lnN          1.10       -          -          -          -          -          -
SigEf_PDF       lnN          1.10       -          -          -          -          -          -
SigEf_PU        lnN          1.10       -          -          -          -          -          -
SigEf_Lep       lnN          1.10       -          -          -          -          -          -
SigEf_Trig      lnN          1.10       -          -          -          -          -          -
NormF_JES       lnN          1.15       -          -          -          -          -          -
NormF_model     lnN          1.15       -          -          -          -          -          -
lumi_13TeV      lnN          1.10       1.10       1.10       1.10       1.10       1.10       1.10
crossSction     lnN          1.10       1.10       1.10       1.10       1.10       1.10       1.10
'''

## template datacard ends here 

## Write templat datacard to the text file with placeholders.
datacard = open('DataCard_MXXXGeV.txt','w')
datacard.write(templatecard)
datacard.close()

## Function to provide the normalization weight factors
def Normalize(n,xs,tot):
    yield_ = n*xs*1./tot
    return yield_

## map of placeholder used in the Template datacard.
## This is analysis specific.
nameinnumber=['ZH',
              'ZZ',
              'WZ',
              'WW',
              'TTBAR',
              'DYJETS',
              'DATA']

## List of signal samples for which limit is needed. 
## This is analysis specific.
signalnameinnumber=[ 'M800',
                    'M1000',
                    'M1200',
                    'M1400',
                    'M1600',
                    'M1800',
                    'M2000',
                    'M2500',
                    'M3000',
                    'M3500',
                    'M4000']

## create the names of place RATE holders
placeholder = [x + "RATE" for x in nameinnumber]
## print placeholder


## valuemap for background and signal with a default value
valuemap = {
    "default" : 0.0
    }

signalvaluemap = {
    "default" : 0.0
    }

## Read the signal background numbers from plain TEXTFile
## this value map is used later to get the datacard by replacing the
## place holders with values stored in this map.
numbers = open(inputtextfilename,'r')
for iline in numbers:
    a,b = iline.split()
    for iname in range(len(nameinnumber)):
        if a==nameinnumber[iname]:
            stringtoprint = nameinnumber[iname]+" value is "+b
            print stringtoprint
            ratename = nameinnumber[iname]+"RATE"
            valuemap[ratename]=b
    ### Following lines fill the 
    ### value map for signal points
    for isigname in range(len(signalnameinnumber)):
        if a==signalnameinnumber[isigname]:
            stringtoprint = signalnameinnumber[isigname]+" value is "+b
            print stringtoprint
            ratename = signalnameinnumber[isigname]+"RATE"
            signalvaluemap[ratename]=b

print valuemap
print signalvaluemap

## Method to access the rootfiles
## Use it to clone and then 
#sigTFile = ROOT.TFile('Merged_DMHistosSpring15_1/main-NCUGlobalTuples_M1500.root','READ')
#sigEvent  = sigTFile.Get('CutFlowAndEachCut/h_cutflow_0')
#sigTEvent = sigTFile.Get('nEvents')
#print sigEvent.GetBinContent(7)
#scaledsig = Normalize(sigEvent.GetBinContent(7), SignalXS['M1500'],sigTEvent.GetEntries())
#print scaledsig

def MakeDataCard(masspoint):
    datacard = open('DataCard_MXXXGeV.txt','r')
    newdatacardname = dirtosave+'/DataCard_'+masspoint+'GeV_MonoHbb_13TeV.txt'
#    os.system('rm '+newdatacardname)
    datacard600 = open(newdatacardname,'w')
    
    for line in datacard:
        ## replace the input root file
        line = line.replace('INPUTROOTFILE', inputrootfilename) 

        ## replace the background values
        for ival in range(len(placeholder)):
            line = line.replace(placeholder[ival],valuemap[placeholder[ival]])
        
        ## replace the signal values
        masspointrate = masspoint + "RATE"
        line = line.replace('SIGNALRATE', signalvaluemap[masspointrate])

        ## replace the signal names
        massname = 'SIG'+masspoint
        line = line.replace('SIGNAL', massname)

        datacard600.write(line)
    datacard600.close()

for imasspoint in range(len(signalnameinnumber)):
    MakeDataCard(signalnameinnumber[imasspoint])

print "datacards produced"
