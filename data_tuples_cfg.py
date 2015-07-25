import FWCore.ParameterSet.Config as cms 

process = cms.Process('ANA')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

#process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
process.GlobalTag.globaltag = 'GR_P_V56'

# HCAL filter (next 2 lines) added on 24 July 2015
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.TFileService=cms.Service("TFileService",
        fileName=cms.string("ntuple_output.root"),
        closeFileFast = cms.untracked.bool(True)
)

process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(True),
)

process.out = cms.OutputModule('PoolOutputModule',
        fileName = cms.untracked.string('edm_output.root'),
        outputCommands = cms.untracked.vstring('keep *')
)

#remove all MC matching when running on data
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'])


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'root://eoscms//eos/cms/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FE26BEB8-D575-E411-A13E-00266CF2AE10.root'
    #'root://eoscms//eos/cms/store/data/Run2015B/JetHT/MINIAOD/PromptReco-v1/000/251/163/00000/F05CF208-A026-E511-85F4-02163E011B15.root',
    #'root://eoscms//eos/cms/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/162/00000/12284DB9-4227-E511-A438-02163E013674.root'
    'root://eoscms//eos/cms/store/data/Run2015B/JetHT/MINIAOD/PromptReco-v1/000/251/163/00000/F05CF208-A026-E511-85F4-02163E011B15.root'
	)
)

#luminosity
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = 'JSON_CRAB/json_DCSONLY_Run2015B.txt').getVLuminosityBlockRange()

#input file
#import FWCore.Utilities.FileUtils as FileUtils
#files2015data = FileUtils.loadListFromFile ('files.txt') 
#readFiles = cms.untracked.vstring( *files2015data )
#process.source.fileNames = readFiles

# Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules_el = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff']

#add them to the VID producer
for idmod in my_id_modules_el:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

switchOnVIDPhotonIdProducer(process, dataFormat)
my_id_modules_ph = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']

#add them to the VID producer
for idmod in my_id_modules_ph:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.bhana = cms.EDAnalyzer('BHAnalyzerTLBSM',
  beamSpot = cms.InputTag('offlineBeamSpot'),
  electronTag = cms.InputTag("slimmedElectrons"),
  muonTag = cms.untracked.InputTag("slimmedMuons"),
  jetTag = cms.untracked.InputTag("slimmedJets"),
  tauTag = cms.untracked.InputTag("slimmedTaus"),
  metTag = cms.untracked.InputTag("slimmedMETs"),
  photonTag = cms.InputTag("slimmedPhotons"),
  rho_lable    = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
  ebRecHitTag = cms.untracked.InputTag("reducedEgamma", "reducedEBRecHits"),
  eeRecHitTag = cms.untracked.InputTag("reducedEgamma", "reducedEERecHits"),
  primaryVertex = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
  triggerTag = cms.untracked.InputTag("TriggerResults","","HLT"),
  
  verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
  conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),

  eleVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
  eleLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
  eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
  eleTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
 
  phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
  phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium"),
  phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight"),
 
  MCLabel = cms.untracked.bool(True),                               
  DEBUG = cms.untracked.bool(False)                               
)


process.p = cms.Path(process.HBHENoiseFilterResultProducer* #produces HBHE bools
	process.ApplyBaselineHBHENoiseFilter*  #reject events based 
	(process.egmPhotonIDSequence+process.egmGsfElectronIDSequence) * process.bhana)

