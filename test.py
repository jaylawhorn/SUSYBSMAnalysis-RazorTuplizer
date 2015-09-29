
from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = "test_AOD_tupler_v4"
config.General.workArea = "crab"

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "razorTuplizer_MC_25ns.py"
config.JobType.allowUndistributedCMSSW = False

config.section_("Data")
config.Data.ignoreLocality = True
config.Data.inputDataset = "/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v3/AODSIM"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.totalUnits = 10
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
config.Data.outLFNDirBase = '/store/user/jlawhorn/'
