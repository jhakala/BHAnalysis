from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'BHnTuples_QCD_HT2000toInf'
config.General.workArea = 'crab_jobs_BHnTuples_QCD_HT2000toInf_Dec7'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'qcdMC.py'
config.JobType.inputFiles=['Summer15_25nsV6_MC.db','Summer15_25nsV6_MC']

config.Data.inputDataset = '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 35
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'BHnTuples_QCD_HT2000toInf_Dec7'

config.Site.storageSite = 'T3_US_Brown'
