from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'BHnTuples_2016G'
config.General.workArea = 'crab_jobs_2016G_Brown'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'bh80Xcfg.py'

config.Data.inputDataset = '/JetHT/Run2016G-PromptReco-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 50
config.Data.lumiMask='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-283059_13TeV_PromptReco_Collisions16_JSON_NoL1T_MuonPhys.txt'
config.Data.outLFNDirBase = '/store/user/%s' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'BHnTuples_2016G'

#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T3_US_Brown'
