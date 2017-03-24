from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'BHnTuples_2015D_76X_May25'
config.General.workArea = 'crab_jobs_2015D_76X_May25'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'bh76cfg_Run2015D_25ns-2May2016.py'
#config.JobType.psetName = 'pmptRecoV4_tuples_2015D_v4.py'
#config.JobType.inputFiles=['Summer15_25nsV6_DATA.db','Summer15_25nsV6_DATA']

config.Data.inputDataset = '/JetHT/Run2015D-16Dec2015-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 35
config.Data.lumiMask='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'BHnTuples_2015D_76x_May25'

config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T3_US_Brown'
