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
runOnData=True #data/MC switch
usePrivateSQlite=False #use external JECs (sqlite file)
useHFCandidates=True #create an additionnal NoHF slimmed MET collection if the option is set to false
applyResiduals=True #application of residual corrections. Have to be set to True once the 13 TeV residual 
#corrections are available. False to be kept meanwhile. Can be kept to False later for private tests or 
#for analysis checks and developments (not the official recommendation!).
#===================================================================
#from Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff import *
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if runOnData:
  process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v4'
else:
  #process.GlobalTag.globaltag = 'auto:run2_mc'
  process.GlobalTag.globaltag = 'MCRUN2_74_V9'

#### For applying jet/met corrections from a sql
#if usePrivateSQlite:
#  from CondCore.DBCommon.CondDBSetup_cfi import *
#  import os
#  if runOnData:
#    jecfile="Summer15_25nsV5_DATA"
#  else:
#    jecfile="Summer15_50nsV4_MC"
# # dBFile can be called by following two ways
# # dBFile = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/PatAlgos/test/"+jecfile+".db")
# # dBFile = os.path.expandvars("/afs/cern.ch/work/a/asaddiqu/BH_CMS/CMSSW_7_4_6_patch1/src/BH_CMS2015/BHAnalysis-master/jec/"+jecfile+".db")
#  #dBFile = os.path.expandvars(jecfile+".db") #A try for crab
#  process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
#    connect = cms.string( "sqlite_file:/afs/cern.ch/user/j/johakala/work/public/CMSSW_7_4_12_patch4/src/BH/BHAnalysis/Summer15_25nsV5_DATA.db" ),
#    toGet =  cms.VPSet(
#    cms.PSet(
#      record = cms.string("JetCorrectionsRecord"),
#      tag = cms.string("JetCorrectorParametersCollection_"+jecfile+"_AK4PF"),
#      label= cms.untracked.string("AK4PF")
#      ),
#    cms.PSet(
#      record = cms.string("JetCorrectionsRecord"),
#      tag = cms.string("JetCorrectorParametersCollection_"+jecfile+"_AK4PFchs"),
#      label= cms.untracked.string("AK4PFchs")
#      ),
#    )
#    )
#process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

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
#-------------------------------------------
#uncertainty file (also can be called by following two ways)
#jecUncertaintyFile="BH/BHAnalysis/Summer15_25nsV5_DATA/Summer15_25nsV5_DATA_UncertaintySources_AK4PFchs.txt"
#jecUncertaintyFile="PhysicsTools/PatUtils/data/Summer15_50nsV4_DATA_UncertaintySources_AK4PFchs.txt"

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

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

# Define the input source
#if runOnData:
#  fname = '

#else:
#  fname = 'root://eoscms.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/60000/001C7571-0511-E511-9B8E-549F35AE4FAF.root'
# Define the input source
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring( 
#'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v4/000/258/177/00000/CE175343-706D-E511-9957-02163E0142B1.root'
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
### the lines below remove the L2L3 residual corrections when processing data
### -------------------------------------------------------------------
if not applyResiduals:
  process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
  process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
  process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
  process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
  process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
  process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
  #if not useHFCandidates:
  #  process.patPFMetT1T2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
  #  process.patPFMetT1T2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
  #  process.patPFMetT2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
  #  process.patPFMetT2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
  #  process.shiftedPatJetEnDownNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
  #  process.shiftedPatJetEnUpNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
### ------------------------------------------------------------------
#luminosity
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
if runOnData:
	process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()

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
my_id_modules_el = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']

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
  filterTag = cms.untracked.InputTag("TriggerResults","","RECO"),
  prescales = cms.InputTag("patTrigger"), 
  verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
  conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),

  eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
  eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
  eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
  eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
 
 # TODO: These need to be updated when Run2 25ns cut-based photonID comes available
  phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
  phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium"),
  phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight"),
 
  MCLabel = cms.untracked.bool(True),                               
  DEBUG = cms.untracked.bool(False)                               
)


process.p = cms.Path(
  process.HBHENoiseFilterResultProducer * # get HBHENoiseFilter decisions
  process.ApplyBaselineHBHENoiseFilter *  # filter based on HBHENoiseFilter decisions
  process.ApplyHBHEIsoNoiseFilter *       # filter for HBHENoise isolation
#  process.eeBadScFilter *                 # apply the EE bad supercrystal filter
  (process.egmPhotonIDSequence+process.egmGsfElectronIDSequence) *
  process.bhana
)
process.p +=cms.Sequence(process.JEC)
