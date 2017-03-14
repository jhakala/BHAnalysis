import os

def AddLines(crab,pSet):
	dataSetName = crab["dataSetName"]
	GlobalTag   = crab["GlobalTag"]

	lines = []
	lines.append("config.General.requestName   = 'BHnTuples_%s'\n"     %dataSetName )
	lines.append("config.General.workArea      = 'crab_jobs_%s'\n"     %dataSetName )
	lines.append("config.Data.inputDataset     = '/JetHT/%s/MINIAOD'\n"%dataSetName )
	lines.append("config.Data.outputDatasetTag = 'BHnTuples_%s'\n"     %dataSetName )
	lines.append("config.JobType.psetName      =  '%s'\n"              % pSet)
	lines.append("config.JobType.pyCfgParams   = ['GlobalTag=%s']\n"   % GlobalTag  )
	return lines

def makeConfigFromTemplate(crab, template, pSet):
	if not os.path.exists(template):
		print "%s does not exits"%template
		return ""
	dataSetName = crab["dataSetName"]
	tempF  = open(template,"r")
	configName = "CrabConfig_%s.py"%dataSetName
	config     = open(configName,"w")
	for line in tempF:
		config.write(line)
	tempF.close()
	for line in AddLines(crab,pSet):
		config.write(line)
	print "	Writing %s"%configName
	return configName


template     = "crabConfig_template.py"
pSet         = "/afs/cern.ch/user/k/kakwok/work/public/Blackhole/CMSSW_8_0_26_patch1/src/BH/BHAnalysis/bh80Xcfg.py"
crabs = []
crabs.append({"dataSetName":"Run2016B-03Feb2017_ver2-v2","GlobalTag":"80X_dataRun2_2016SeptRepro_v7"})
crabs.append({"dataSetName":"Run2016C-03Feb2017-v1"     ,"GlobalTag":"80X_dataRun2_2016SeptRepro_v7"})
crabs.append({"dataSetName":"Run2016D-03Feb2017-v1"     ,"GlobalTag":"80X_dataRun2_2016SeptRepro_v7"})
crabs.append({"dataSetName":"Run2016E-03Feb2017-v1"     ,"GlobalTag":"80X_dataRun2_2016SeptRepro_v7"})
crabs.append({"dataSetName":"Run2016F-03Feb2017-v1"     ,"GlobalTag":"80X_dataRun2_2016SeptRepro_v7"})
crabs.append({"dataSetName":"Run2016G-03Feb2017-v1"     ,"GlobalTag":"80X_dataRun2_2016SeptRepro_v7"})
crabs.append({"dataSetName":"Run2016H-03Feb2017_ver2-v1","GlobalTag":"80X_dataRun2_Prompt_v16"})
crabs.append({"dataSetName":"Run2016H-03Feb2017_ver3-v1","GlobalTag":"80X_dataRun2_Prompt_v16"})

for crab in crabs:
	CrabConfig = makeConfigFromTemplate(crab,template,pSet)
	print "Submitting %s" % CrabConfig
	os.system("crab submit %s" % CrabConfig)
