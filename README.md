# BHAnalysis
##1)CMSSW VERSION
-------------
CMSSW_7_4_5 version adopted becausue photon/electron IDs are implemented with this version.
##2) Required packages 
-----------------
for photon/electron IDs, do the followings:
```
cmsrel CMSSW_7_4_5
cd CMSSW_7_4_5/src
cmsenv
scram b -j 10  
```
##3) Reference Twikis for IDS
---------------------------
Jet ID     : https://twiki.cern.ch/twiki/bin/view/CMS/TopJME  
Electron ID: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2  
Photon ID  : https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2  
Muon ID    : https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO  
##4) Cuts Applied in Ntuplizer
-----------------------
###Cuts for all Objects:
if (|eta|>2.4 || pt < 20 GeV) Reject these jets/electrons/photons/muons.

###Additional cuts for Electron:
if(!e->gsfTrack()) Reject these electrons

###Additional cuts for Muon:
if(!mu->isGlobalMuon() || !mu->isPFMuon())Reject these muons
