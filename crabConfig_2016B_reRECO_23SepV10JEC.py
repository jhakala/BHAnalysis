from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'BHnTuples_2016B_reRECO_23SepV10JEC'
config.General.workArea = 'crab_jobs_2016B_reRECO_23SepV10JEC'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'bh80Xcfg.py'
config.JobType.inputFiles=['Spring16_23Sep2016AllV1_DATA.db']

config.Data.inputDataset = '/JetHT/Run2016B-23Sep2016-v3/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 50
config.Data.lumiMask='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.outLFNDirBase = '/store/user/%s' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'BHnTuples_2016B_reRECO_23SepV10JEC'

#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T3_US_Brown'
