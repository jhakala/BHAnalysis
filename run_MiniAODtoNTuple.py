import os
import glob
#'/afs/cern.ch/user/k/kakwok/eos/cms/store/user/kakwok/MiniAOD/Charybdis_BH10_CH_MD4000_MBH9000_n6/Charybdis_BH10_CH_MD4000_MBH9000_n6/160817_010102/0000/miniAOD-prod_PAT_1.root'
Topdir="/afs/cern.ch/user/k/kakwok/eos/cms/store/user/kakwok/MiniAOD/"
#Outdir="/mnt/hadoop/users/mkwok/NTuple/"
folder_list = os.listdir(Topdir)
print folder_list


jobdate="1608"
for dir in folder_list:
	#print "%s%s/*/%s*/*/*.root"%(Topdir,dir,jobdate)
	flist = glob.glob("%s%s/*/%s*/*/*.root"%(Topdir,dir,jobdate))
	fname = flist[0]
	#outname = Outdir+dir+"_NTuple.root"
	outname = dir+"_NTuple.root"
	cmd = "cmsRun bh80Xcfg_MC.py InputFile=%s OutputFile=%s"%( fname, outname)
	if not(os.path.exists(outname)):
	        print cmd
	#        os.system(cmd)
	else:
	        print "Already produced this masspoint %s...skipping" % dir
