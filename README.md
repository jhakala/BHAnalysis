# BHAnalysis
##1) Instructions
-------------------------------------------
###i) Compile the nTuplizer against CMSSW_8_0_26_patch1 (or later)
```
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
cmsenv
git cms-init
git cms-merge-topic -u cms-met:METRecipe_8020
git cms-merge-topic ikrav:egm_id_80X_v3
scram b -j12
mkdir BH
cd BH
git clone https://github.com/kakwok/BHAnalysis.git
scram b -j8
```
Updated for ReMiniAOD campaign.
###iia) Customize bh80Xcfg.py 
For a particular task, the relevant things to customize are:
- runOnData -- should be True for data, False for MC.
- globaltag -- make sure you apply the correct global tag, see Useful Links below
- lumisToProcess -- for running on data, this typically points to a "golden JSON" from the data certificiation team.
- process.source -- for running locally (as opposed to on the GRID), this is the input file for your task.
- fileName -- for running locally, this is the output file for your task.
- process.maxevents -- for running locally, this is the number of events to run over. (-1 for all events)
- eleLooseIdMap, phoMediumIdMap, etc. -- these should use the most up-to-date ID maps with the correct bunch crossing spacing.

note: As this python script needs to be updated often, look for the most updated ```*cfg.py```

###iib) Customize crabConfig_2016G.py
The relevant things to customize for running on the grid are:
- inputdataSet -- this should point to the dataset as found on DAS and MiniAOD campaign twiki, see Useful Links below.
- lumiMask -- for processing data, this should point toward the golden JSON.
- storageSite -- this should point to your T2 or T3 where you have write access.

###iiia) Run on data/mc locally
```
cd BHAnalysis
cmsenv
cmsRun bh80Xcfg.py GlobalTag=80X_dataRun2_2016SeptRepro_v7
```
Check which global tag to use from Global Tags and MiniAOD
To run a list of signal samples, 
```
vim run_MiniAODtoNTuple.py
python run_MiniAODtoNTyple.py
```
This uses bh80Xcfg_MC.py

###iiib) Run on data/mc over the grid using CRAB3

```
cd BHAnalysis
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crabConfig_2016G.py
```
##2) Useful links:
* Data Aggregation Service : https://cmsweb.cern.ch/das/
* Global Tags              : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
* JetMet                   : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription  
* Jet ID                   : https://twiki.cern.ch/twiki/bin/view/CMS/TopJME  
* Electron ID              : https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2  
* Photon ID                : https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2  
* Muon ID                  : https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO  
* MET Filters              : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#MiniAOD
* HBHEnoiseFilters         : https://twiki.cern.ch/twiki/bin/viewauth/CMS/HCALNoiseFilterRecipe
* GoldenJSON               : https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/
* MiniAOD                  : https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
