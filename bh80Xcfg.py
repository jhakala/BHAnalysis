import FWCore.ParameterSet.Config as cms 
import FWCore.ParameterSet.VarParsing as VarParsing
#from RecoMET.METFilters.eeBadScFilter_cfi import *

process = cms.Process('ANA')

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.Geometry_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1

#------------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------------

options = VarParsing.VarParsing()
options.register('GlobalTag',
                "auto",
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "GlobalTag to be used")
options.register('outputFile',
                "ntuple_output.root",
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "filename of output root file")

options.parseArguments()

print "Going to use GlobalTag = %s"% options.GlobalTag

#------------------------------------------------------------------------------------
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

#===MET filters in 80X ==============================================
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
#For tagging mode, i.e. saving the decision
process.BadPFMuonFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
#For tagging mode, i.e. saving the decision
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

#process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')
#======================================================================

# Bad EE supercrystal filter
#process.load(eeBadScFilter)


#configurable options ==============================================
runOnData=True #data/MC switch
usePrivateSQlite=False #use external JECs (sqlite file)
useHFCandidates=True #create an additionnal NoHF slimmed MET collection if the option is set to false
applyResiduals=True #application of residual JES corrections. Setting this to false removes the residual JES corrections.
#===================================================================

#==Global tags ====================================================
#see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
#===================================================================

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if runOnData:
  #process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v14'           # For Prompt-RECO 2016H only
  #process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7' # For reMiniAOD
  process.GlobalTag.globaltag = options.GlobalTag 
else:
  #process.GlobalTag.globaltag = 'auto:run2_mc'
  process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_RealisticBS_25ns_13TeV2016_v1_mc'
#===================================================================

#==For applying Jet energy correction from a sqlite file ======================
# from: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecSqliteFile
# 1. Get sqlite file from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
# 2. Find list of tags by  conddb --db<dbfile.db> listTags
# 3. Update the db name
#===================================================================
if usePrivateSQlite:
  process.load("CondCore.DBCommon.CondDBCommon_cfi")
  from CondCore.DBCommon.CondDBSetup_cfi import *
  process.jec = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
      messageLevel = cms.untracked.int32(0)
      ),
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(
    cms.PSet(
        record = cms.string('JetCorrectionsRecord'),
        tag    = cms.string('JetCorrectorParametersCollection_Spring16_23Sep2016AllV2_DATA_AK4PF'),
        label  = cms.untracked.string('AK4PF')
        ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag    = cms.string('JetCorrectorParametersCollection_Spring16_23Sep2016AllV2_DATA_AK4PFchs'),
      label  = cms.untracked.string('AK4PFchs')
      ),
    ), 
    connect = cms.string('sqlite:Spring16_23Sep2016AllV2_DATA.db')
  )
  process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

#==Update JEC after MINIAOD ====================================================================
# from: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets
#==============================================================================================
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

if runOnData:
    updateJetCollection(
       process,
       jetSource = cms.InputTag('slimmedJets'),
       labelName = 'UpdatedJEC',
       jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']), 'None')  # Do not forget 'L2L3Residual' on data!
    )
else:
    updateJetCollection(
       process,
       jetSource = cms.InputTag('slimmedJets'),
       labelName = 'UpdatedJEC',
       jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')  # Do not forget 'L2L3Residual' on data!
    )
#==============================================================================================
# Configure output option
#==============================================================================================

process.load('FWCore.MessageService.MessageLogger_cfi')

process.TFileService=cms.Service("TFileService",
        fileName=cms.string(options.outputFile),
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

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring( 
#'file:/afs/cern.ch/user/k/kakwok/eos/cms/store/data/Run2016B/JetHT/MINIAOD/23Sep2016-v3/00000/00144F9E-BA97-E611-A8B0-00259074AE48.root'
#'file:/afs/cern.ch/user/k/kakwok/work/public/CMSSW_7_6_5/src/Blackhole/BHAnalysis/eos/cms/store/data/Run2015C_25ns/JetHT/MINIAOD/16Dec2015-v1/20000/D41FEE23-49B5-E511-B288-3417EBE6471D.root'
#'file:/afs/cern.ch/user/k/kakwok/eos/cms/store/data/Run2016C/JetHT/MINIAOD/PromptReco-v2/000/275/890/00000/B08F2A69-5A3F-E611-BA56-02163E01477C.root'
#'file:/afs/cern.ch/user/k/kakwok/eos/cms/store/data/Run2016H/JetHT/MINIAOD/PromptReco-v2/000/281/256/00000/CEF4A29D-6E82-E611-8CF7-02163E01215C.root'
'file:/afs/cern.ch/user/k/kakwok/eos/cms/store/data/Run2016B/JetHT/MINIAOD/03Feb2017_ver2-v2/110000/003A92CA-6FED-E611-82CD-0025905B8590.root'
#'file:/afs/cern.ch/user/k/kakwok/eos/cms/store/data/Run2016B/JetHT/MINIAOD/23Sep2016-v3/00000/00144F9E-BA97-E611-A8B0-00259074AE48.root'
#'file:/afs/cern.ch/user/k/kakwok/work/public/Blackhole/CMSSW_8_1_0_pre16/src/BH/BHAnalysis/BH2016G_badEvents_MINIAOD_reRECO.root'
 )
)
# How many events to process
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

#==============================================================================================


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
#if runOnData:
#	process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt').getVLuminosityBlockRange()

# Set up electron ID (VID framework)
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules_el = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']

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
  muonTag = cms.InputTag("slimmedMuons"),
  #jetTag =  cms.InputTag("slimmedJets"),
  jetTag =  cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
  tauTag =  cms.InputTag("slimmedTaus"),
  metTag =  cms.InputTag("slimmedMETs"),
  photonTag = cms.InputTag("slimmedPhotons"),
  rho_lable    = cms.InputTag("fixedGridRhoFastjetAll"),
  ebRecHitTag = cms.untracked.InputTag("reducedEgamma", "reducedEBRecHits"),
  eeRecHitTag = cms.untracked.InputTag("reducedEgamma", "reducedEERecHits"),
  primaryVertex = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
  badChHadfilter = cms.InputTag("BadChargedCandidateFilter"),
  badMufilter    = cms.InputTag("BadPFMuonFilter"),
  triggerTag = cms.InputTag("TriggerResults","","HLT"),
  filterTag = cms.InputTag("TriggerResults","","PAT"),
  prescales = cms.InputTag("patTrigger"), 
  verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
  conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),

  eleVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"  ),
  eleLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose" ),
  eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
  eleTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight" ),
 
  phoLooseIdMap  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose" ),
  phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
  phoTightIdMap  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight" ),
 
  MCLabel = cms.untracked.bool(False),                               
  DEBUG = cms.untracked.bool(False)                               
)


process.p = cms.Path(
  process.HBHENoiseFilterResultProducer * # get HBHENoiseFilter decisions
  process.ApplyBaselineHBHENoiseFilter *  # filter based on HBHENoiseFilter decisions
  process.ApplyHBHEIsoNoiseFilter *       # filter for HBHENoise isolation
  (process.egmPhotonIDSequence+process.egmGsfElectronIDSequence) *
  process.BadPFMuonFilter *		  # 80x new met filter
  process.BadChargedCandidateFilter *     # 80x new met filter
  process.bhana
)
#process.p +=cms.Sequence(process.JEC)
