# BHAnalysis
##1) INSTRUCTION TO RUN NTUPLIZER ON LXPLUS
-------------------------------------------
###i) Compile nTuplizer against CMSSW_7_4_12_patch4
For running on lxplus, do:
```
cmsrel CMSSW_7_4_12_patch4
mkdir CMSSW_7_4_12_patch4/src/<some dir>
cd CMSSW_7_4_12_patch4/src/<some dir>
git clone https://github.com/jhakala/BHAnalysis.git
scram b -j8
```
###ii) Change path inside maketuples_jec_cfg.py
Please change all the "/afs/cern.ch/work/<your username>" paths according to your lxplus account.
###v) Run for data/mc on lxplus locally
Inside "<some dir>/BHAnalysis/maketuples_jec_cfg.py"   
keep "runOnData=True" for data  
and  "runOnData=False" for MC
and then simply
```
cd BHAnalysis
cmsRun maketuples_jec_cfg.py
```
If you get any error at this point then please let me know.

##2) Useful Twikis:
JetMet     : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription  
Jet ID     : https://twiki.cern.ch/twiki/bin/view/CMS/TopJME  
Electron ID: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2  
Photon ID  : https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2  
Muon ID    : https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO  
