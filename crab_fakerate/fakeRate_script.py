import os
import sys
import optparse
import ROOT
import re

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleRECOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleIDSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonIDISOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.FakeRateProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import *
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles, runsAndLumis
### main python file to run ###

def main():

  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('--year', dest='year', help='which year sample', default='2018', type='string')
  parser.add_option('-m', dest='ismc', help='to apply sf correction or not', default=True, action='store_true')
  parser.add_option('-d', dest='ismc', help='to apply sf correction or not', action='store_false')
  (opt, args) = parser.parse_args()

  if opt.ismc:
    if opt.year == "2016a":
      p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),puWeight_2016_preAPV(),PrefCorr2016(),muonIDISOSF2016apv(),muonScaleRes2016a(),eleRECOSF2016apv(),eleIDSF2016apv(), jmeCorrections_UL2016APVMC(),FakeRate2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016b":
      p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),puWeight_2016_postAPV(),PrefCorr2016(),muonIDISOSF2016(),muonScaleRes2016b(),eleRECOSF2016(),eleIDSF2016(),jmeCorrections_UL2016MC(),FakeRate2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017":
      p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),puWeight_2017(),PrefCorr2017(),muonIDISOSF2017(),muonScaleRes2017(),eleRECOSF2017(),eleIDSF2017(), jmeCorrections_UL2017MC(), FakeRate2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2018":
      p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),puWeight_2018(),muonIDISOSF2018(),muonScaleRes2018(),eleRECOSF2018(),eleIDSF2018(),jmeCorrections_UL2018MC(), FakeRate2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")


# Sequence for data
  if not (opt.ismc):
    if opt.year == "2016b":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016a(),jmeCorrections_UL2016B(),FakeRate2016apv()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016c":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016a(),jmeCorrections_UL2016C(),FakeRate2016apv()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016d":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016a(),jmeCorrections_UL2016D(),FakeRate2016apv()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016e":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016a(),jmeCorrections_UL2016E(),FakeRate2016apv()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016f_apv":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016a(),jmeCorrections_UL2016APVF(),FakeRate2016apv()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016f":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016b(),jmeCorrections_UL2016F(),FakeRate2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016g":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016b(),jmeCorrections_UL2016G(),FakeRate2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016h":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016b(),jmeCorrections_UL2016H(),FakeRate2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017b":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),jmeCorrections_UL2017B(),FakeRate2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017c":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),jmeCorrections_UL2017C(),FakeRate2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017d":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),jmeCorrections_UL2017D(),FakeRate2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017e":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),jmeCorrections_UL2017E(),FakeRate2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017f":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),jmeCorrections_UL2017F(),FakeRate2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2018a":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2018(),jmeCorrections_UL2018A(),FakeRate2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2018b":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2018(),jmeCorrections_UL2018B(),FakeRate2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2018c":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2018(),jmeCorrections_UL2018C(),FakeRate2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2018d":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2018(),jmeCorrections_UL2018D(),FakeRate2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
  p.run()

if __name__ == "__main__":
    sys.exit(main())
