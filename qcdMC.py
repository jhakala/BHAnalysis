import FWCore.ParameterSet.Config as cms 
#from RecoMET.METFilters.eeBadScFilter_cfi import *

process = cms.Process('ANA')

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

# HBHE noise filter (next 2 lines) added on 24 July 2015
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)


process.ApplyHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
    reverseDecision = cms.bool(False)
)

# Bad EE supercrystal filter
#process.load(eeBadScFilter)

# How many events to process
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#configurable options ==============================================
runOnData=False #data/MC switch
usePrivateSQlite=True #use external JECs (sqlite file)
useHFCandidates=True #create an additionnal NoHF slimmed MET collection if the option is set to false
applyResiduals=True #application of residual JES corrections. Setting this to false removes the residual JES corrections.
#===================================================================
#from Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff import *
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if runOnData:
  process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v4'
else:
  #process.GlobalTag.globaltag = 'auto:run2_mc'
  process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_realisticBS_v1'

#### For applying jet/met corrections from a sql
if usePrivateSQlite:
  process.load("CondCore.DBCommon.CondDBCommon_cfi")
  from CondCore.DBCommon.CondDBSetup_cfi import *
  if runOnData:
    tagname = 'JetCorrectorParametersCollection_Summer15_25nsV6_DATA_AK4PF'
  else:
    tagname = 'JetCorrectorParametersCollection_Summer15_25nsV6_MC_AK4PF'
  process.jec = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
      messageLevel = cms.untracked.int32(0)
      ),
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(
    cms.PSet(
        record = cms.string('JetCorrectionsRecord'),
        tag    = cms.string(tagname),
        label  = cms.untracked.string('AK4PF')
        ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag    = cms.string(tagname+'chs'),
      label  = cms.untracked.string('AK4PFchs')
      ),
    ), 
    connect = cms.string('sqlite:Summer15_25nsV6_MC.db')
  )
  process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

#-----------------For JEC-----------------
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.patJetCorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
  src = cms.InputTag("slimmedJets"),
  levels = ['L1FastJet', 'L2Relative', 'L3Absolute'],
  payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = process.patJetsUpdated.clone(
  jetSource = cms.InputTag("slimmedJets"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
  )
process.JEC = cms.Sequence( process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC )

process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",
        fileName=cms.string("ntuple_output.root"),
        closeFileFast = cms.untracked.bool(True)
)

process.out = cms.OutputModule('PoolOutputModule',
  compressionLevel = cms.untracked.int32(4),
  compressionAlgorithm = cms.untracked.string('LZMA'),
  eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
  fileName = cms.untracked.string('edm_output.root'),
  outputCommands = cms.untracked.vstring('keep *'),
  dataset = cms.untracked.PSet(
          filterName = cms.untracked.string(''),
          dataTier = cms.untracked.string('')
          ),
  dropMetaData = cms.untracked.string('ALL'),
  fastCloning = cms.untracked.bool(False),
  overrideInputFileSplitLevels = cms.untracked.bool(True)
)

if runOnData:
  process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring( 
  #'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/177/00000/CE175343-706D-E511-9957-02163E0142B1.root'
   )
  )
else:
  process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring( 
  'root://eoscms.cern.ch//store/mc/RunIISpring15DR74/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/10000/64929C34-C018-E511-8FB4-002481E0D5CE.root'
   )
  )


### ---------------------------------------------------------------------------
### Removing the HF from the MET computation
### ---------------------------------------------------------------------------
#if not useHFCandidates:
#  process.noHFCands = cms.EDFilter("CandPtrSelector",
#  src=cms.InputTag("packedPFCandidates"),
#  cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
#  )
#jets are rebuilt from those candidates by the tools, no need to do anything else
### =================================================================================
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

#default configuration for miniAOD reprocessing, change the isData flag to run on data
#for a full met computation, remove the pfCandColl input
#runMetCorAndUncFromMiniAOD(process,
#  isData=runOnData,
  #jecUncFile=jecUncertaintyFile
#  )

#if not useHFCandidates:
#  runMetCorAndUncFromMiniAOD(process,
#  isData=runOnData,
#  pfCandColl=cms.InputTag("noHFCands"),
#  jecUncFile=jecUncertaintyFile,
#  postfix="NoHF"
#  )
### -------------------------------------------------------------------
### remove the L2L3 residual corrections when processing data
### -------------------------------------------------------------------
#if not applyResiduals:
#  process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
#  process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
#  process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
#  process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
#  process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
#  process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
  #if not useHFCandidates:
  #  process.patPFMetT1T2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
  #  process.patPFMetT1T2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
  #  process.patPFMetT2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
  #  process.patPFMetT2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
  #  process.shiftedPatJetEnDownNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
  #  process.shiftedPatJetEnUpNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
### ------------------------------------------------------------------
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
if runOnData:
  process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()

# Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules_el = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']

#add them to the VID producer
for idmod in my_id_modules_el:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

switchOnVIDPhotonIdProducer(process, dataFormat)
my_id_modules_ph = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff']

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
  filterTag = cms.untracked.InputTag("TriggerResults","","RECO"),
  prescales = cms.InputTag("patTrigger"), 
  verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
  conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),

  eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
  eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
  eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
  eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
 
  phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
  phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
  phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
 
  MCLabel = cms.untracked.bool(not runOnData),                               
  DEBUG = cms.untracked.bool(False)                               
)


if runOnData:
  process.p = cms.Path(
    process.HBHENoiseFilterResultProducer * # get HBHENoiseFilter decisions
    process.ApplyBaselineHBHENoiseFilter *  # filter based on HBHENoiseFilter decisions
    process.ApplyHBHEIsoNoiseFilter *       # filter for HBHENoise isolation
    (process.egmPhotonIDSequence+process.egmGsfElectronIDSequence) *
    process.bhana
  )
else:
  process.p = cms.Path(
    (process.egmPhotonIDSequence+process.egmGsfElectronIDSequence) *
    process.bhana
  )
process.p +=cms.Sequence(process.JEC)
