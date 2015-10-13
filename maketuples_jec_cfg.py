import FWCore.ParameterSet.Config as cms 

process = cms.Process('ANA')

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

# HCAL filter (next 2 lines) added on 24 July 2015
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)
# How many events to process
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#configurable options ==============================================
runOnData=True #data/MC switch
usePrivateSQlite=True #use external JECs (sqlite file)
useHFCandidates=True #create an additionnal NoHF slimmed MET collection if the option is set to false
applyResiduals=True #application of residual corrections. Have to be set to True once the 13 TeV residual 
#corrections are available. False to be kept meanwhile. Can be kept to False later for private tests or 
#for analysis checks and developments (not the official recommendation!).
#===================================================================
#from Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff import *
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if runOnData:
  process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v2'
else:
  #process.GlobalTag.globaltag = 'auto:run2_mc'
  process.GlobalTag.globaltag = 'MCRUN2_74_V9'
if usePrivateSQlite:
  from CondCore.DBCommon.CondDBSetup_cfi import *
  import os
  if runOnData:
    jecfile="Summer15_25nsV5_DATA"
  else:
    jecfile="Summer15_50nsV4_MC"
 # dBFile can be called by following two ways
 # dBFile = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/PatAlgos/test/"+jecfile+".db")
 # dBFile = os.path.expandvars("/afs/cern.ch/work/a/asaddiqu/BH_CMS/CMSSW_7_4_6_patch1/src/BH_CMS2015/BHAnalysis-master/jec/"+jecfile+".db")
  #dBFile = os.path.expandvars(jecfile+".db") #A try for crab
  process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
    connect = cms.string( "sqlite_file:/afs/cern.ch/user/j/johakala/work/public/CMSSW_7_4_12_patch4/src/BH/BHAnalysis/Summer15_25nsV5_DATA.db" ),
    toGet =  cms.VPSet(
    cms.PSet(
      record = cms.string("JetCorrectionsRecord"),
      tag = cms.string("JetCorrectorParametersCollection_"+jecfile+"_AK4PF"),
      label= cms.untracked.string("AK4PF")
      ),
    cms.PSet(
      record = cms.string("JetCorrectionsRecord"),
      tag = cms.string("JetCorrectorParametersCollection_"+jecfile+"_AK4PFchs"),
      label= cms.untracked.string("AK4PFchs")
      ),
    )
    )
process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')
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
jecUncertaintyFile="BH/BHAnalysis/Summer15_25nsV5_DATA/Summer15_25nsV5_DATA_UncertaintySources_AK4PFchs.txt"
#jecUncertaintyFile="PhysicsTools/PatUtils/data/Summer15_50nsV4_DATA_UncertaintySources_AK4PFchs.txt"

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.TFileService=cms.Service("TFileService",
        fileName=cms.string("ntuple_output_13.root"),
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
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/584/00000/8A5B7CB0-855D-E511-BFCF-02163E0133A9.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/587/00000/F664AC07-935D-E511-A019-02163E01424B.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/629/00000/6EDCF302-0A5F-E511-A5EF-02163E014642.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/630/00000/86ACFECD-3C5F-E511-B8F2-02163E014374.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/662/00000/300565D4-F55E-E511-95AF-02163E011D25.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/672/00000/0C762538-075F-E511-ACEE-02163E0136C1.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/673/00000/F20A98EE-1C5F-E511-8845-02163E014767.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/674/00000/DA9B86E1-F95E-E511-B75E-02163E013460.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/675/00000/6467F0A5-A75F-E511-AE8A-02163E013389.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/675/00000/7285B0A0-A75F-E511-8C73-02163E011911.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/675/00000/D68BBF9F-A75F-E511-9D37-02163E011BD2.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/1069ABF5-C15F-E511-898B-02163E014160.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/206851F5-C15F-E511-A5A5-02163E0142BA.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/2A98DEF6-C15F-E511-B5F9-02163E012A2E.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/36C1F29D-C25F-E511-B94E-02163E011F85.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/3A2126E4-C15F-E511-A4D2-02163E011888.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/506D1AE8-C15F-E511-8A08-02163E0144EA.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/5CC9C1EC-C15F-E511-9417-02163E014370.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/62B3323F-C35F-E511-B49D-02163E0143A2.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/866DBCDD-C15F-E511-B4DE-02163E013965.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/C0DE1AED-C15F-E511-971C-02163E014337.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/C4318ED7-C15F-E511-91C0-02163E012175.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/EC93F0E5-C15F-E511-940F-02163E01350C.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/676/00000/EEE0E0DD-C15F-E511-A1FD-02163E014268.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/677/00000/58B82BAD-985F-E511-AA44-02163E012386.root', 
'root://eoscms.cern.ch//store/data/Run2015D/JetHT/MINIAOD/PromptReco-v3/000/256/677/00000/5E2C8E5D-985F-E511-A44D-02163E013394.root' 
 )
)

### ---------------------------------------------------------------------------
### Removing the HF from the MET computation
### ---------------------------------------------------------------------------
if not useHFCandidates:
  process.noHFCands = cms.EDFilter("CandPtrSelector",
  src=cms.InputTag("packedPFCandidates"),
  cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
  )
#jets are rebuilt from those candidates by the tools, no need to do anything else
### =================================================================================
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

#default configuration for miniAOD reprocessing, change the isData flag to run on data
#for a full met computation, remove the pfCandColl input
runMetCorAndUncFromMiniAOD(process,
  isData=runOnData,
  jecUncFile=jecUncertaintyFile
  )

if not useHFCandidates:
  runMetCorAndUncFromMiniAOD(process,
  isData=runOnData,
  pfCandColl=cms.InputTag("noHFCands"),
  jecUncFile=jecUncertaintyFile,
  postfix="NoHF"
  )
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
  if not useHFCandidates:
    process.patPFMetT1T2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT1T2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.shiftedPatJetEnDownNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    process.shiftedPatJetEnUpNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
### ------------------------------------------------------------------
#luminosity
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
if runOnData:
	process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-257599_13TeV_PromptReco_Collisions15_25ns_JSON.txt').getVLuminosityBlockRange()

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
  metTag = cms.untracked.InputTag("slimmedMETsNoHF"),
  photonTag = cms.InputTag("slimmedPhotons"),
  rho_lable    = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
  ebRecHitTag = cms.untracked.InputTag("reducedEgamma", "reducedEBRecHits"),
  eeRecHitTag = cms.untracked.InputTag("reducedEgamma", "reducedEERecHits"),
  primaryVertex = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
  triggerTag = cms.untracked.InputTag("TriggerResults","","HLT"),
  prescales = cms.InputTag("patTrigger"), 
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
process.p +=cms.Sequence(process.JEC)