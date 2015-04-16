import FWCore.ParameterSet.Config as cms 

process = cms.Process('ANA')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process.GlobalTag.globaltag = 'PHYS14_25_V1::All'

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

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

process.out = cms.OutputModule('PoolOutputModule',                                                                                                           fileName = cms.untracked.string('edm_output.root'),                                                                              
        outputCommands = cms.untracked.vstring('keep *')
)

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/mc/Phys14DR/QCD_HT-500To1000_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/827D2760-916F-E411-8F4C-AC162DABAF78.root')
    fileNames = cms.untracked.vstring('root://eoscms//eos/cms/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/FE26BEB8-D575-E411-A13E-00266CF2AE10.root')
)

########### Photon ID ###########
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")

########### Elec Id #############
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
# overwrite a default parameter: for miniAOD, the collection name is a slimmed one
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')

from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)

my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V1_miniAOD_cff']
#Add them to the VID producer
for idmod in my_id_modules:
      setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.bhana = cms.EDAnalyzer('BHAnalyzerTLBSM',
  electronTag = cms.untracked.InputTag("slimmedElectrons"),
  muonTag = cms.untracked.InputTag("slimmedMuons"),
  jetTag = cms.untracked.InputTag("slimmedJets"),
  tauTag = cms.untracked.InputTag("slimmedTaus"),
  metTag = cms.untracked.InputTag("slimmedMETs"),
  photonTag = cms.untracked.InputTag("slimmedPhotons"),
  rho_lable    = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
  ebRecHitTag = cms.untracked.InputTag("reducedEgamma", "reducedEBRecHits"),
  eeRecHitTag = cms.untracked.InputTag("reducedEgamma", "reducedEERecHits"),
  primaryVertex = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
  triggerTag = cms.untracked.InputTag("TriggerResults"),
  full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
  phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
  phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
  phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
  electronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-veto"),
  electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-medium"),
  MCLabel = cms.untracked.bool(True),                               
  DEBUG = cms.untracked.bool(False)                               
)



#process.p = cms.Path(process.bhana)
#process.p = cms.Path(process.photonIDValueMapProducer*process.bhana)
#process.p = cms.Path(process.egmGsfElectronIDSequence * process.bhana)
process.p = cms.Path((process.photonIDValueMapProducer+ process.egmGsfElectronIDSequence) * process.bhana)

