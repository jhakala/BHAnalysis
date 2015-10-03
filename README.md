# BHAnalysis
##1) INSTRUCTION TO RUN NTUPLIZER ON LXPLUS
-------------------------------------------
###i) Put Ntuplizer in src directory
Assuming you have a directory "BHAnalysis" somewhere in your lxplus account containing Ntuplizer, which is cloned from git. Then please follow this in src directory
```
mkdir BH_CMS2015
cp -r /afs/cern.ch/...etc../../BHAnalysis BH_CMS2015/
scramv1 b -j 8
cd BH_CMS2015/BHAnalysis
```
###ii) Change path inside maketuples_jec_cfg.py
Please change all the "/afs/cern.ch/work/a/asaddiqu/" paths according to your lxplus account.
###v) Run for data/mc on lxplus locally
Inside "maketuples_jec_cfg.py"   
keep "runOnData=True" for data  
and  "runOnData=False" for MC
and then simply
```
cmsRun maketuples_jec_cfg.py
```
If you get any error at this point then please let me know.

##2) Useful Twikis:
JetMet     : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription  
Jet ID     : https://twiki.cern.ch/twiki/bin/view/CMS/TopJME  
Electron ID: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2  
Photon ID  : https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2  
Muon ID    : https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO  
