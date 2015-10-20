# BHAnalysis
##1) INSTRUCTION TO RUN NTUPLIZER ON LXPLUS
-------------------------------------------
###i) Compile the nTuplizer against CMSSW_7_4_14
For running on lxplus, do:
```
cmsrel CMSSW_7_4_14
mkdir CMSSW_7_4_14/src/<some dir>
cd CMSSW_7_4_14/src/<some dir>
git clone https://github.com/jhakala/BHAnalysis.git
scram b -j8
```
###iia) Customize pmptRecoV4_tuples_2015D.py
For a particular task, the relevant things to customize are:
  > runOnData -- should be True for data, False for MC.
  > globaltag -- make sure you apply the correct global tag, see Useful Links below
  > lumisToProcess -- for running on data, this typically points to a "golden JSON" from the data certificiation team.
  > process.source -- for running locally (as opposed to on the GRID), this is the input file for your task.
  > fileName -- for running locally, this is the output file for your task.
  > process.maxevents -- for running locally, this is the number of events to run over. (-1 for all events)
  > eleLooseIdMap, phoMediumIdMap, etc. -- these should use the most up-to-date ID maps with the correct bunch crossing spacing.

###iib) Customize crabConfig_2015DpmptRecoV4.py
The relevant things to customize for running on the grid are:
  > inputdataSet -- this should point to the dataset as found on DAS, see Useful Links below.
  > lumiMask -- for processing data, this should point toward the golden JSON.
  > storageSite -- this should point to your T2 or T3 where you have write access.

###iiia) Run for data/mc on lxplus locally
```
cd BHAnalysis
cmsenv
cmsRun pmptRecoV4_tuples_2015D.py
```

###iiib) Run for data/mc on lxplus on the grid using CRAB3

```
cd BHAnalysis
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crabConfig_2015DpmptRecoV4.py
```
##2) Useful links:
Data Aggregation Service : https://cmsweb.cern.ch/das/
Global Tags              : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
JetMet                   : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription  
Jet ID                   : https://twiki.cern.ch/twiki/bin/view/CMS/TopJME  
Electron ID              : https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2  
Photon ID                : https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2  
Muon ID                  : https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO  
MET Filters              : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#MiniAOD
