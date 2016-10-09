import sys
import os
import csv
#import ROOT
#ROOT.gROOT.SetBatch(True)  
#from DataSetInfo import *

if len(sys.argv) < 6 :
    print "insufficient options provided see help function "
    exit (1)

if len(sys.argv) == 6 :
    print ('You are making datacards for '+sys.argv[1]+' with '+sys.argv[2]+' and datacards will be saved in '+sys.argv[5])

inputtextfilename=sys.argv[1]
inputrootfilename=sys.argv[2]
inputsignaluncfile=sys.argv[3]
inputotheruncfile=sys.argv[4]
dirtosave=sys.argv[5]

os.system('mkdir -p '+dirtosave)

## prepare template datacard and store in this variable. 
## do whatever change you want to do 
## and keep variables in capital letters with prefix "T" e.g. signal strength can be written as TSIGNAL and DYJets can be TDYJETS
## these can be replaced by either a python or shell script to make datacards for a specific mass point and analysis. 
templatecard='''
imax    1        number of channels
jmax    *        number of backgrounds
kmax    *        number of nuisance parameters (sources of systematical uncertainties)

---------------------------------------------------------------------------

shapes  *        ZPTOZHMM  INPUTROOTFILE  $PROCESS  $PROCESS_$SYSTEMATIC

---------------------------------------------------------------------------

bin              ZPTOZHMM
observation      DATARATE

---------------------------------------------------------------------------

bin                          ZPTOZHMM     ZPTOZHMM     ZPTOZHMM  
process                      SIGNAL       ZJETS        SUBDOM
process                      0            1            2         
rate                         SIGNALRATE   ZJETSRATE    SUBDOMRATE

---------------------------------------------------------------------------

background_bTag     shape    -            1.000        -              
background_QCD      shape    -            1.000        -              
background_PDF      shape    -            1.000        -              
background_JES      shape    -            1.000        -              
background_FitDev   shape    -            1.000        - 
background_FitGood  lnN      -            FITGOOD      -      
background_Norm     lnN      -            NORM         -
background_SubDom   lnN      -            -            SUBDOM
signal_bTag         lnN      SIGBTAG      -            -              
signal_QCD          lnN      SIGQCD       -            -              
signal_PDF          lnN      SIGPDF       -            -              
signal_JES          lnN      SIGJES       -            -              
signal_PU           lnN      SIGPU        -            -              
signal_Lep          lnN      SIGLEP       -            -              
signal_Trig         lnN      SIGTRIG      -            -              
lumi_13TeV          lnN      1.027        1.027        1.027
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
nameinnumber=['SUBDOM',
              'ZJETS',
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

## Read the signal background numbers from plain TEXT File
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

# Read the uncertainty numbers according mass point and sources
def sigUncValue(source,masspoint):
    myData = csv.DictReader(open(inputsignaluncfile), delimiter="\t")
    for row in myData:
        if row['mass'] == masspoint:
            return row[source]

# Read the uncertainty numbers 
def otherUncValue(source):
    myData = csv.DictReader(open(inputotheruncfile), delimiter=" ")
    for row in myData:
        return row[source]

# Make datacards
def MakeDataCard(masspoint):
    datacard = open('DataCard_MXXXGeV.txt','r')
    newdatacardname = dirtosave+'/DataCard_'+masspoint+'GeV_MonoHbb_13TeV.txt'
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
        
        mass = masspoint.replace('M','')
        ## replace the uncertainty values
        line = line.replace('SIGBTAG', sigUncValue('bTag',mass))
        line = line.replace('SIGQCD',  sigUncValue('QCD',mass))
        line = line.replace('SIGPDF',  sigUncValue('PDF',mass))
        line = line.replace('SIGJES',  sigUncValue('JES',mass))
        line = line.replace('SIGPU',   sigUncValue('pileUp',mass))
        line = line.replace('SIGLEP',  sigUncValue('lepton',mass))
        line = line.replace('SIGTRIG', sigUncValue('Trigger',mass))

        line = line.replace('FITGOOD', otherUncValue('fitGood'))
        line = line.replace('NORM',    otherUncValue('norm'))
        line = line.replace('SUBDOM',  otherUncValue('subDom'))

        datacard600.write(line)
    datacard600.close()

for imasspoint in range(len(signalnameinnumber)):
    MakeDataCard(signalnameinnumber[imasspoint])
    
print "datacards produced"
