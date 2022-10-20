from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = 'QCD_Muenriched_120to170'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
# hadd nano will not be needed once nano tools are in cmssw
config.JobType.inputFiles = ['fakeRate_script.py', '../scripts/haddnano.py','keep_and_drop.txt']
config.JobType.sendPythonFolder = True

config.section_("Data")
config.Data.inputDataset = '/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL17NanoAODv2-106X_mc2017_realistic_v8-v1/NANOAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.outLFNDirBase = '/store/group/phys_top/ExtraYukawa/'
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.outputDatasetTag = 'QCD_Muenriched_120to170'

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
