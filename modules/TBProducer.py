import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import numpy as np
from numpy import sign

class TBProducer(Module):
  def __init__(self , year):
    self.year = year
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("HLT_passEle32WPTight", "I")
    self.out.branch("lhe_nlepton", "I")
    self.out.branch("n_tight_muon", "I")
    self.out.branch("n_fakeable_muon", "I")
    self.out.branch("n_loose_muon", "I")
    self.out.branch("n_tight_ele", "I")
    self.out.branch("n_fakeable_ele", "I")
    self.out.branch("n_loose_ele", "I")
    self.out.branch("n_tight_jet", "I")
    self.out.branch("n_bjet_DeepB_M", "I")
    self.out.branch("n_bjet_DeepB_L", "I")
    self.out.branch("n_tight_nob", "I")
    self.out.branch("HT", "F")
    self.out.branch("nHad_tau", "I")
    self.out.branch("j1_pt", "F")
    self.out.branch("j1_eta", "F")
    self.out.branch("j1_phi", "F")
    self.out.branch("j1_mass", "F")
    self.out.branch("j1_isB", "F")
    self.out.branch("j2_pt", "F")
    self.out.branch("j2_eta", "F")
    self.out.branch("j2_phi", "F")
    self.out.branch("j2_mass", "F")
    self.out.branch("j2_isB", "F")
    self.out.branch("j3_pt", "F")
    self.out.branch("j3_eta", "F")
    self.out.branch("j3_phi", "F")
    self.out.branch("j3_mass", "F")
    self.out.branch("j3_isB", "F")
    self.out.branch("j4_pt", "F")
    self.out.branch("j4_eta", "F")
    self.out.branch("j4_phi", "F")
    self.out.branch("j4_mass", "F")
    self.out.branch("j4_isB", "F")
    self.out.branch("mj1j2", "F")
    self.out.branch("mj1j3", "F")
    self.out.branch("mj1j4", "F")
    self.out.branch("mj2j3", "F")
    self.out.branch("mj2j4", "F")
    self.out.branch("mj3j4", "F")
    self.out.branch("mj1j2j3", "F")
    self.out.branch("mj1j2j4", "F")
    self.out.branch("mj2j3j4", "F")
    self.out.branch("mj1j2j3j4", "F")
    self.out.branch("drj1j2", "F")
    self.out.branch("drj1j3", "F")
    self.out.branch("drj1j4", "F")
    self.out.branch("drj2j3", "F")
    self.out.branch("drj2j4", "F")
    self.out.branch("drj3j4", "F")

    self.out.branch("j1_nob_pt", "F")
    self.out.branch("j1_nob_eta", "F")
    self.out.branch("j1_nob_phi", "F")
    self.out.branch("j1_nob_mass", "F")
    self.out.branch("j2_nob_pt", "F")
    self.out.branch("j2_nob_eta", "F")
    self.out.branch("j2_nob_phi", "F")
    self.out.branch("j2_nob_mass", "F")
    self.out.branch("j3_nob_pt", "F")
    self.out.branch("j3_nob_eta", "F")
    self.out.branch("j3_nob_phi", "F")
    self.out.branch("j3_nob_mass", "F")
    self.out.branch("j4_nob_pt", "F")
    self.out.branch("j4_nob_eta", "F")
    self.out.branch("j4_nob_phi", "F")
    self.out.branch("j4_nob_mass", "F")
    self.out.branch("mj1j2_nob", "F")

    self.out.branch("j1_b_pt", "F")
    self.out.branch("j1_b_eta", "F")
    self.out.branch("j1_b_phi", "F")
    self.out.branch("j1_b_mass", "F")
    self.out.branch("j2_b_pt", "F")
    self.out.branch("j2_b_eta", "F")
    self.out.branch("j2_b_phi", "F")
    self.out.branch("j2_b_mass", "F")
    self.out.branch("j3_b_pt", "F")
    self.out.branch("j3_b_eta", "F")
    self.out.branch("j3_b_phi", "F")
    self.out.branch("j3_b_mass", "F")
    self.out.branch("j4_b_pt", "F")
    self.out.branch("j4_b_eta", "F")
    self.out.branch("j4_b_phi", "F")
    self.out.branch("j4_b_mass", "F")
    self.out.branch("mj1j2_b", "F")

    self.out.branch("Onelep_region", "B")
    self.out.branch("Twolep_region", "B")
    self.out.branch("fourJet_region", "B")

    self.out.branch("Onelep_l1_faketag", "B")
    self.out.branch("Onelep_l1_pt", "F")
    self.out.branch("Onelep_l1_eta", "F")
    self.out.branch("Onelep_l1_phi", "F")
    self.out.branch("Onelep_l1_mass", "F")
    self.out.branch("Onelep_l1_id", "I")
    self.out.branch("Onelep_l1_pdgid", "I")
    self.out.branch("Onelep_l1_isprompt", "I")

    self.out.branch("Twolep_channel", "I")
    self.out.branch("Twolep_2P0F", "B")
    self.out.branch("Twolep_1P1F", "B")
    self.out.branch("Twolep_0P2F", "B")
    self.out.branch("Twolep_l1_faketag", "B")
    self.out.branch("Twolep_l1_pt", "F")
    self.out.branch("Twolep_l1_eta", "F")
    self.out.branch("Twolep_l1_phi", "F")
    self.out.branch("Twolep_l1_mass", "F")
    self.out.branch("Twolep_l1_id", "I")
    self.out.branch("Twolep_l1_pdgid", "I")
    self.out.branch("Twolep_l1_isprompt", "I")
    self.out.branch("Twolep_l2_faketag", "B")
    self.out.branch("Twolep_l2_pt", "F")
    self.out.branch("Twolep_l2_eta", "F")
    self.out.branch("Twolep_l2_phi", "F")
    self.out.branch("Twolep_l2_mass", "F")
    self.out.branch("Twolep_l2_id", "I")
    self.out.branch("Twolep_l2_pdgid", "I")
    self.out.branch("Twolep_l2_isprompt", "I")
    self.out.branch("Twolep_mll", "F")

    self.out.branch("tightJets_id","I",lenVar="nJet")
    self.out.branch("tightJets_nob_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_DeepCSVmedium_id","I",lenVar="nJet")
    self.out.branch("tightJets_b_DeepCSVloose_id","I",lenVar="nJet")
    self.out.branch("tightElectrons_id","I",lenVar="nElectron")
    self.out.branch("fakeable_Electrons_id","I",lenVar="nElectron")
    self.out.branch("additional_vetoElectrons_id","I",lenVar="nElectron")
    self.out.branch("tightMuons_id","I",lenVar="nMuon")
    self.out.branch("fakeable_Muons_id","I",lenVar="nMuon")
    self.out.branch("additional_looseMuons_id","I",lenVar="nMuon")
    self.out.branch("Had_tau_id","I",lenVar="nTau")

    self.out.branch("tb_region","B")
    self.out.branch("tt1L_region","B")
    self.out.branch("wjet_region","B")
    self.out.branch("wv_region","B")
    self.out.branch("tt2L_region","B")
    self.out.branch("dy_region","B")
    self.out.branch("met_user","F")
    self.out.branch("met_phi_user","F")

    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    self.is_lhe = bool(inputTree.GetBranch("nLHEPart"))

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    # PV selection
    if (event.PV_npvsGood<1): return False

    # trigger selection
    # special action for 2017 single ele HLT, https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations#Single_Electron_Triggers
    HLT_passEle32WPTight=0
    if self.year=="2017":
      trgobjs=Collection(event, 'TrigObj')
      if event.HLT_Ele32_WPTight_Gsf_L1DoubleEG==1:
	for iobj in range(0,event.nTrigObj):
	  if trgobjs[iobj].id==11 and (trgobjs[iobj].filterBits & (1<<10))== (1<<10):
	    HLT_passEle32WPTight=1

    self.out.fillBranch("HLT_passEle32WPTight",HLT_passEle32WPTight)

    lhe_nlepton=0
    if self.is_lhe:
      lheparticle = Collection(event, 'LHEPart')
      for ilhe in range(0, event.nLHEPart):
        if lheparticle[ilhe].status==1 and (abs(lheparticle[ilhe].pdgId)==11 or abs(lheparticle[ilhe].pdgId)==13 or abs(lheparticle[ilhe].pdgId)==15):
          lhe_nlepton=lhe_nlepton+1

    self.out.fillBranch("lhe_nlepton", lhe_nlepton)

    # total number of ele+muon, currently require at least 1 leptons
    if ((event.nMuon + event.nElectron) < 1): return False
    if not event.nJet>2: return False

    # Muon selection: tight cut-based ID + tight PF iso, or loose cut-based ID + loose PF iso, with pt > 20 GeV
    muons = Collection(event, 'Muon')
    muon_v4_temp=TLorentzVector()
    muon_v4_temp_raw=TLorentzVector()
    tightMuons = []
    tightMuons_raw = []
    tightMuons_pdgid = []
    tightMuons_id = []
    fakeable_Muons = []
    fakeable_Muons_pdgid = []
    fakeable_Muons_id = []
    additional_looseMuons = []
    additional_looseMuons_pdgid = []
    additional_looseMuons_id = []

    jet_v4_temp=TLorentzVector()

    for imu in range(0, event.nMuon):
      # following cuts are preseletion for MVA muon ID
      if abs(muons[imu].eta)>2.4 or muons[imu].pfRelIso04_all>0.4: continue
      if not (muons[imu].looseId and event.Muon_corrected_pt[imu]>15):continue
      
      if muons[imu].tightId:
        if muons[imu].pfRelIso04_all<0.15:
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          muon_v4_temp_raw.SetPtEtaPhiM(muons[imu].pt, muons[imu].eta, muons[imu].phi, muons[imu].mass)
          tightMuons.append(muon_v4_temp.Clone())
          tightMuons_raw.append(muon_v4_temp_raw.Clone())
          tightMuons_pdgid.append(muons[imu].pdgId)
          tightMuons_id.append(imu)
        elif muons[imu].pfRelIso04_all>0.2:
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          fakeable_Muons.append(muon_v4_temp.Clone())
          fakeable_Muons_pdgid.append(muons[imu].pdgId)
          fakeable_Muons_id.append(imu)
        else:
          muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
          additional_looseMuons.append(muon_v4_temp.Clone())
          additional_looseMuons_pdgid.append(muons[imu].pdgId)
          additional_looseMuons_id.append(imu)
      else:
        muon_v4_temp.SetPtEtaPhiM(event.Muon_corrected_pt[imu], muons[imu].eta, muons[imu].phi, muons[imu].mass)
        additional_looseMuons.append(muon_v4_temp.Clone())
        additional_looseMuons_pdgid.append(muons[imu].pdgId)
        additional_looseMuons_id.append(imu)

    n_tight_muon = len(tightMuons)
    n_fakeable_muon = len(fakeable_Muons)
    n_loose_muon = len(additional_looseMuons)

    self.out.fillBranch("n_tight_muon", n_tight_muon)
    self.out.fillBranch("n_fakeable_muon", n_fakeable_muon)
    self.out.fillBranch("n_loose_muon", n_loose_muon)
    tightMuons_id.extend(np.zeros(event.nMuon-len(tightMuons_id),int)-1)
    fakeable_Muons_id.extend(np.zeros(event.nMuon-len(fakeable_Muons_id),int)-1)
    additional_looseMuons_id.extend(np.zeros(event.nMuon-len(additional_looseMuons_id),int)-1)
    self.out.fillBranch("tightMuons_id", tightMuons_id)
    self.out.fillBranch("fakeable_Muons_id", fakeable_Muons_id)
    self.out.fillBranch("additional_looseMuons_id", additional_looseMuons_id)

    # electron selection: tight (veto) cut-based ID + impact parameter cut, with pt > 15 GeV
    electrons = Collection(event, 'Electron')
    electron_v4_temp=TLorentzVector()
    electron_v4_temp_raw=TLorentzVector()
    tightElectrons = []
    tightElectrons_raw = []
    tightElectrons_pdgid = []
    tightElectrons_id = []
    fakeable_Electrons = []
    fakeable_Electrons_pdgid = []
    fakeable_Electrons_id = []
    additional_vetoElectrons = []
    additional_vetoElectrons_pdgid = []
    additional_vetoElectrons_id = []


    for iele in range(0, event.nElectron):
      # following cuts are preseletion for MVA electron ID
#      if abs(electrons[iele].eta)>2.5 or abs(electrons[iele].dxy)>0.05 or abs(electrons[iele].dz)>0.1: continue
      if not (abs(electrons[iele].eta)<2.5 and electrons[iele].cutBased>0 and electrons[iele].pt>15):continue

      if (abs(electrons[iele].deltaEtaSC+electrons[iele].eta)<1.479 and abs(electrons[iele].dxy)<0.05 and abs(electrons[iele].dz)<0.1) or (abs(electrons[iele].deltaEtaSC+electrons[iele].eta)>1.479 and abs(electrons[iele].dxy)<0.1 and abs(electrons[iele].dz)<0.2):
        if electrons[iele].cutBased>3:
          electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          electron_v4_temp_raw.SetPtEtaPhiM(electrons[iele].pt/electrons[iele].eCorr, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass/electrons[iele].eCorr)
          tightElectrons.append(electron_v4_temp.Clone())
          tightElectrons_raw.append(electron_v4_temp_raw.Clone())
          tightElectrons_pdgid.append(electrons[iele].pdgId)
          tightElectrons_id.append(iele)
        elif electrons[iele].cutBased>1:
          electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          fakeable_Electrons.append(electron_v4_temp.Clone())
          fakeable_Electrons_pdgid.append(electrons[iele].pdgId)
          fakeable_Electrons_id.append(iele)
        else:
          electron_v4_temp.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
          additional_vetoElectrons.append(electron_v4_temp.Clone())
          additional_vetoElectrons_pdgid.append(electrons[iele].pdgId)
          additional_vetoElectrons_id.append(iele)


    n_tight_ele = len(tightElectrons)
    n_fakeable_ele = len(fakeable_Electrons)
    n_loose_ele = len(additional_vetoElectrons)
    self.out.fillBranch("n_tight_ele", n_tight_ele)
    self.out.fillBranch("n_fakeable_ele", n_fakeable_ele)
    self.out.fillBranch("n_loose_ele", n_loose_ele)
    tightElectrons_id.extend(np.zeros(event.nElectron-len(tightElectrons_id),int)-1)
    fakeable_Electrons_id.extend(np.zeros(event.nElectron-len(fakeable_Electrons_id),int)-1)
    additional_vetoElectrons_id.extend(np.zeros(event.nElectron-len(additional_vetoElectrons_id),int)-1)
    self.out.fillBranch("tightElectrons_id", tightElectrons_id)
    self.out.fillBranch("fakeable_Electrons_id", fakeable_Electrons_id)
    self.out.fillBranch("additional_vetoElectrons_id", additional_vetoElectrons_id)

    # tight leptons and additional loose leptons collection
    tightLeptons = tightMuons + tightElectrons
    tightLeptons.sort(key=lambda x: x.Pt(), reverse=True)
    tightLeptons_raw = tightMuons_raw + tightElectrons_raw
    tightLeptons_raw.sort(key=lambda x: x.Pt(), reverse=True)
    fakeableLeptons = fakeable_Muons + fakeable_Electrons
    fakeableLeptons.sort(key=lambda x: x.Pt(), reverse=True)
    looseLeptons = additional_looseMuons + additional_vetoElectrons
    looseLeptons.sort(key=lambda x: x.Pt(), reverse=True)

    tau_v4_temp=TLorentzVector()
    taus = Collection(event, 'Tau')
    nHad_tau=0
    Had_tau_id=[]
    for itau in range(0, event.nTau):
      tau_v4_temp.SetPtEtaPhiM(taus[itau].pt, taus[itau].eta, taus[itau].phi, taus[itau].mass)
      pass_tau_lep_Dr=1
      if taus[itau].pt>20 and abs(taus[itau].eta)<2.3 and taus[itau].idDecayModeOldDMs and taus[itau].idDeepTau2017v2p1VSe>=4 and taus[itau].idDeepTau2017v2p1VSjet>=4 and taus[itau].idDeepTau2017v2p1VSmu>=1:
	for ilep in range(0,len(tightLeptons)):
          if tau_v4_temp.DeltaR(tightLeptons[ilep])<0.4:pass_tau_lep_Dr=0
	for ilep in range(0,len(fakeableLeptons)):
          if tau_v4_temp.DeltaR(fakeableLeptons[ilep])<0.4:pass_tau_lep_Dr=0
	if pass_tau_lep_Dr:
	  nHad_tau=nHad_tau+1
	  Had_tau_id.append(itau)
    self.out.fillBranch("nHad_tau", nHad_tau)

    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
    # tightLepVeto PF jets (ak4), 2016 (111=7), 2017/2018 (110=6), medium B-tag WP
    # DeepCSV=(nanoaod btagDeepB) loose: 0.1355, medium: 0.4506, tight: 0.7738
    # DeepFlavor=(nanoaod btagDeepFlavB) loose: 0.0532, medium: 0.3040, tight: 0.7476

    # c-jet tag is based on two-D cuts, medium DeepJet WP:
    # CvsL=btagDeepFlavCvL: 0.085, CvsB=btagDeepFlavCvB: 0.34
    # c-tag not available in NANOAOD yet

    jets = Collection(event, 'Jet')

    j1_pt=-99
    j1_eta=-99
    j1_phi=-99
    j1_mass=-99
    j1_isB=-99
    j2_pt=-99
    j2_eta=-99
    j2_phi=-99
    j2_mass=-99
    j2_isB=-99
    j3_pt=-99
    j3_eta=-99
    j3_phi=-99
    j3_mass=-99
    j3_isB=-99
    j4_pt=-99
    j4_eta=-99
    j4_phi=-99
    j4_mass=-99
    j4_isB=-99
    mj1j2=-99
    mj1j3=-99
    mj1j4=-99
    mj2j3=-99
    mj2j4=-99
    mj3j4=-99
    mj1j2j3=-99
    mj1j2j4=-99
    mj2j3j4=-99
    mj1j2j3j4=-99
    drj1j2=-99
    drj1j3=-99
    drj1j4=-99
    drj2j3=-99
    drj2j4=-99
    drj3j4=-99

    tightJets_id = []
    tightJets_nob_id = []

    tightJets_b_DeepCSVmedium_id = []
    tightJets_b_DeepCSVloose_id = []

    # require DeltaR between Jets and tight leptons greater than 0.4
    jet_v4_all = []
    bjet_v4_all = []
    nobjet_v4_all = []
    for ijet in range(0, event.nJet):

      if abs(jets[ijet].eta)>4.7:continue

      jet_is_tau=0
      if nHad_tau>0:
        for ita in Had_tau_id:
          if ijet==event.Tau_jetIdx[ita]:jet_is_tau=1
      if jet_is_tau:continue

      pass_jet_lep_Dr=1
      jet_v4_temp.SetPtEtaPhiM(event.Jet_pt_nom[ijet],event.Jet_eta[ijet],event.Jet_phi[ijet],event.Jet_mass_nom[ijet])
      for ilep in range(0,len(tightLeptons)):
	if jet_v4_temp.DeltaR(tightLeptons[ilep])<0.4:pass_jet_lep_Dr=0
      for ilep in range(0,len(fakeableLeptons)):
	if jet_v4_temp.DeltaR(fakeableLeptons[ilep])<0.4:pass_jet_lep_Dr=0
      for ilep in range(0,len(looseLeptons)):
	if jet_v4_temp.DeltaR(looseLeptons[ilep])<0.4:pass_jet_lep_Dr=0

      if not (pass_jet_lep_Dr>0):continue
      if not (jets[ijet].jetId==6 and event.Jet_pt_nom[ijet]>30):continue 

      tightJets_id.append(ijet)
      jet_v4_all.append(jet_v4_temp.Clone())

      if abs(jets[ijet].eta)<2.5:

        if self.year=="2016apv":
          if jets[ijet].btagDeepFlavB > 0.2598:
            tightJets_b_DeepCSVmedium_id.append(ijet)
            bjet_v4_all.append(jet_v4_temp.Clone())
          elif jets[ijet].btagDeepFlavB > 0.0508:
            tightJets_b_DeepCSVloose_id.append(ijet)
            bjet_v4_all.append(jet_v4_temp.Clone())
  
        if self.year=="2016":
          if jets[ijet].btagDeepFlavB > 0.2489:
            tightJets_b_DeepCSVmedium_id.append(ijet)
            bjet_v4_all.append(jet_v4_temp.Clone())
          elif jets[ijet].btagDeepFlavB > 0.0480:
            tightJets_b_DeepCSVloose_id.append(ijet)
            bjet_v4_all.append(jet_v4_temp.Clone())
  
        if self.year=="2017":
          if jets[ijet].btagDeepFlavB > 0.3040:
            tightJets_b_DeepCSVmedium_id.append(ijet)
            bjet_v4_all.append(jet_v4_temp.Clone())
          elif jets[ijet].btagDeepFlavB > 0.0532:
            tightJets_b_DeepCSVmedium_id.append(ijet)
            bjet_v4_all.append(jet_v4_temp.Clone())
  
        if self.year=="2018":
          if jets[ijet].btagDeepFlavB > 0.2783:
            tightJets_b_DeepCSVmedium_id.append(ijet)
            bjet_v4_all.append(jet_v4_temp.Clone())
          elif jets[ijet].btagDeepFlavB > 0.0490:
            tightJets_b_DeepCSVloose_id.append(ijet)
            bjet_v4_all.append(jet_v4_temp.Clone())

    HT=0
    for ijet in range(0,len(tightJets_id)):
      HT=HT+event.Jet_pt_nom[tightJets_id[ijet]]
    self.out.fillBranch("HT",HT)

    tightJets_nob_id = [x for x in tightJets_id if (x not in tightJets_b_DeepCSVmedium_id and x not in tightJets_b_DeepCSVloose_id)]
    nobjet_v4_all = [x for x in jet_v4_all if x not in bjet_v4_all]

    n_tight_jet = len(tightJets_id)
    n_bjet_DeepB_M = len(tightJets_b_DeepCSVmedium_id)
    n_bjet_DeepB_L = len(tightJets_b_DeepCSVloose_id)
    n_tight_nob = len(tightJets_nob_id)
    self.out.fillBranch("n_tight_jet",n_tight_jet)
    self.out.fillBranch("n_bjet_DeepB_M",n_bjet_DeepB_M)
    self.out.fillBranch("n_bjet_DeepB_L",n_bjet_DeepB_L)
    self.out.fillBranch("n_tight_nob",n_tight_nob)

    Had_tau_id.extend(np.zeros(event.nTau-len(Had_tau_id),int)-1)
    self.out.fillBranch("Had_tau_id", Had_tau_id)
    

    if n_tight_jet>3:
      j4_pt=event.Jet_pt_nom[tightJets_id[3]]
      j4_eta=event.Jet_eta[tightJets_id[3]]
      j4_phi=event.Jet_phi[tightJets_id[3]]
      j4_mass=event.Jet_mass_nom[tightJets_id[3]]
      j4_isB = 1 if (tightJets_id[3] in tightJets_b_DeepCSVmedium_id or tightJets_id[3] in tightJets_b_DeepCSVloose_id) else 0
      j3_pt=event.Jet_pt_nom[tightJets_id[2]]
      j3_eta=event.Jet_eta[tightJets_id[2]]
      j3_phi=event.Jet_phi[tightJets_id[2]]
      j3_mass=event.Jet_mass_nom[tightJets_id[2]]
      j3_isB = 1 if (tightJets_id[2] in tightJets_b_DeepCSVmedium_id or tightJets_id[2] in tightJets_b_DeepCSVloose_id) else 0
      j2_pt=event.Jet_pt_nom[tightJets_id[1]]
      j2_eta=event.Jet_eta[tightJets_id[1]]
      j2_phi=event.Jet_phi[tightJets_id[1]]
      j2_mass=event.Jet_mass_nom[tightJets_id[1]]
      j2_isB = 1 if (tightJets_id[1] in tightJets_b_DeepCSVmedium_id or tightJets_id[1] in tightJets_b_DeepCSVloose_id) else 0
      j1_pt=event.Jet_pt_nom[tightJets_id[0]]
      j1_eta=event.Jet_eta[tightJets_id[0]]
      j1_phi=event.Jet_phi[tightJets_id[0]]
      j1_mass=event.Jet_mass_nom[tightJets_id[0]]
      j1_isB = 1 if (tightJets_id[0] in tightJets_b_DeepCSVmedium_id or tightJets_id[0] in tightJets_b_DeepCSVloose_id) else 0
      mj1j2=(jet_v4_all[0]+jet_v4_all[1]).M()
      mj1j3=(jet_v4_all[0]+jet_v4_all[2]).M()
      mj1j4=(jet_v4_all[0]+jet_v4_all[3]).M()
      mj2j3=(jet_v4_all[1]+jet_v4_all[2]).M()
      mj2j4=(jet_v4_all[1]+jet_v4_all[3]).M()
      mj3j4=(jet_v4_all[2]+jet_v4_all[3]).M()
      mj1j2j3=(jet_v4_all[0]+jet_v4_all[1]+jet_v4_all[2]).M()
      mj1j2j4=(jet_v4_all[0]+jet_v4_all[1]+jet_v4_all[3]).M()
      mj2j3j4=(jet_v4_all[1]+jet_v4_all[2]+jet_v4_all[3]).M()
      mj1j2j3j4=(jet_v4_all[0]+jet_v4_all[1]+jet_v4_all[2]+jet_v4_all[3]).M()
      drj1j2=jet_v4_all[0].DeltaR(jet_v4_all[1])
      drj1j3=jet_v4_all[0].DeltaR(jet_v4_all[2])
      drj1j4=jet_v4_all[0].DeltaR(jet_v4_all[3])
      drj2j3=jet_v4_all[1].DeltaR(jet_v4_all[2])
      drj2j4=jet_v4_all[1].DeltaR(jet_v4_all[3])
      drj3j4=jet_v4_all[2].DeltaR(jet_v4_all[3])

    if n_tight_jet==3:
      j3_pt=event.Jet_pt_nom[tightJets_id[2]]
      j3_eta=event.Jet_eta[tightJets_id[2]]
      j3_phi=event.Jet_phi[tightJets_id[2]]
      j3_mass=event.Jet_mass_nom[tightJets_id[2]]
      j3_isB = 1 if (tightJets_id[2] in tightJets_b_DeepCSVmedium_id or tightJets_id[2] in tightJets_b_DeepCSVloose_id) else 0
      j2_pt=event.Jet_pt_nom[tightJets_id[1]]
      j2_eta=event.Jet_eta[tightJets_id[1]]
      j2_phi=event.Jet_phi[tightJets_id[1]]
      j2_mass=event.Jet_mass_nom[tightJets_id[1]]
      j2_isB = 1 if (tightJets_id[1] in tightJets_b_DeepCSVmedium_id or tightJets_id[1] in tightJets_b_DeepCSVloose_id) else 0
      j1_pt=event.Jet_pt_nom[tightJets_id[0]]
      j1_eta=event.Jet_eta[tightJets_id[0]]
      j1_phi=event.Jet_phi[tightJets_id[0]]
      j1_mass=event.Jet_mass_nom[tightJets_id[0]]
      j1_isB = 1 if (tightJets_id[0] in tightJets_b_DeepCSVmedium_id or tightJets_id[0] in tightJets_b_DeepCSVloose_id) else 0
      mj1j2=(jet_v4_all[0]+jet_v4_all[1]).M()
      mj1j3=(jet_v4_all[0]+jet_v4_all[2]).M()
      mj2j3=(jet_v4_all[1]+jet_v4_all[2]).M()
      mj1j2j3=(jet_v4_all[0]+jet_v4_all[1]+jet_v4_all[2]).M()
      drj1j2=jet_v4_all[0].DeltaR(jet_v4_all[1])
      drj1j3=jet_v4_all[0].DeltaR(jet_v4_all[2])
      drj2j3=jet_v4_all[1].DeltaR(jet_v4_all[2])
    if n_tight_jet==2:
      j2_pt=event.Jet_pt_nom[tightJets_id[1]]
      j2_eta=event.Jet_eta[tightJets_id[1]]
      j2_phi=event.Jet_phi[tightJets_id[1]]
      j2_mass=event.Jet_mass_nom[tightJets_id[1]]
      j2_isB = 1 if (tightJets_id[1] in tightJets_b_DeepCSVmedium_id or tightJets_id[1] in tightJets_b_DeepCSVloose_id) else 0
      j1_pt=event.Jet_pt_nom[tightJets_id[0]]
      j1_eta=event.Jet_eta[tightJets_id[0]]
      j1_phi=event.Jet_phi[tightJets_id[0]]
      j1_mass=event.Jet_mass_nom[tightJets_id[0]]
      j1_isB = 1 if (tightJets_id[0] in tightJets_b_DeepCSVmedium_id or tightJets_id[0] in tightJets_b_DeepCSVloose_id) else 0
      mj1j2=(jet_v4_all[0]+jet_v4_all[1]).M()
      drj1j2=jet_v4_all[0].DeltaR(jet_v4_all[1])
    if n_tight_jet==1:
      j1_pt=event.Jet_pt_nom[tightJets_id[0]]
      j1_eta=event.Jet_eta[tightJets_id[0]]
      j1_phi=event.Jet_phi[tightJets_id[0]]
      j1_mass=event.Jet_mass_nom[tightJets_id[0]]
      j1_isB = 1 if (tightJets_id[0] in tightJets_b_DeepCSVmedium_id or tightJets_id[0] in tightJets_b_DeepCSVloose_id) else 0


    self.out.fillBranch("j1_pt",j1_pt)
    self.out.fillBranch("j1_eta",j1_eta)
    self.out.fillBranch("j1_phi",j1_phi)
    self.out.fillBranch("j1_mass",j1_mass)
    self.out.fillBranch("j1_isB",j1_isB)
    self.out.fillBranch("j2_pt",j2_pt)
    self.out.fillBranch("j2_eta",j2_eta)
    self.out.fillBranch("j2_phi",j2_phi)
    self.out.fillBranch("j2_mass",j2_mass)
    self.out.fillBranch("j2_isB",j2_isB)
    self.out.fillBranch("j3_pt",j3_pt)
    self.out.fillBranch("j3_eta",j3_eta)
    self.out.fillBranch("j3_phi",j3_phi)
    self.out.fillBranch("j3_mass",j3_mass)
    self.out.fillBranch("j3_isB",j3_isB)
    self.out.fillBranch("j4_pt",j4_pt)
    self.out.fillBranch("j4_eta",j4_eta)
    self.out.fillBranch("j4_phi",j4_phi)
    self.out.fillBranch("j4_mass",j4_mass)
    self.out.fillBranch("j4_isB",j4_isB)
    self.out.fillBranch("mj1j2", mj1j2)
    self.out.fillBranch("mj1j3", mj1j3)
    self.out.fillBranch("mj1j4", mj1j4)
    self.out.fillBranch("mj2j3", mj2j3)
    self.out.fillBranch("mj2j4", mj2j4)
    self.out.fillBranch("mj3j4", mj3j4)
    self.out.fillBranch("mj1j2j3", mj1j2j3)
    self.out.fillBranch("mj1j2j4", mj1j2j4)
    self.out.fillBranch("mj2j3j4", mj2j3j4)
    self.out.fillBranch("mj1j2j3j4", mj1j2j3j4)
    self.out.fillBranch("drj1j2", drj1j2)
    self.out.fillBranch("drj1j3", drj1j3)
    self.out.fillBranch("drj1j4", drj1j4)
    self.out.fillBranch("drj2j3", drj2j3)
    self.out.fillBranch("drj2j4", drj2j4)
    self.out.fillBranch("drj3j4", drj3j4)

    j1_nob_pt=-99
    j1_nob_eta=-99
    j1_nob_phi=-99
    j1_nob_mass=-99
    j2_nob_pt=-99
    j2_nob_eta=-99
    j2_nob_phi=-99
    j2_nob_mass=-99
    j3_nob_pt=-99
    j3_nob_eta=-99
    j3_nob_phi=-99
    j3_nob_mass=-99
    j4_nob_pt=-99
    j4_nob_eta=-99
    j4_nob_phi=-99
    j4_nob_mass=-99
    mj1j2_nob=-99

    j1_b_pt=-99
    j1_b_eta=-99
    j1_b_phi=-99
    j1_b_mass=-99
    j2_b_pt=-99
    j2_b_eta=-99
    j2_b_phi=-99
    j2_b_mass=-99
    j3_b_pt=-99
    j3_b_eta=-99
    j3_b_phi=-99
    j3_b_mass=-99
    j4_b_pt=-99
    j4_b_eta=-99
    j4_b_phi=-99
    j4_b_mass=-99
    mj1j2_b=-99

    if (n_bjet_DeepB_M+n_bjet_DeepB_L)>3:
      j1_b_pt=bjet_v4_all[0].Pt()
      j1_b_eta=bjet_v4_all[0].Eta()
      j1_b_phi=bjet_v4_all[0].Phi()
      j1_b_mass=bjet_v4_all[0].M()
      j2_b_pt=bjet_v4_all[1].Pt()
      j2_b_eta=bjet_v4_all[1].Eta()
      j2_b_phi=bjet_v4_all[1].Phi()
      j3_b_mass=bjet_v4_all[2].M()
      j3_b_pt=bjet_v4_all[2].Pt()
      j3_b_eta=bjet_v4_all[2].Eta()
      j3_b_phi=bjet_v4_all[2].Phi()
      j4_b_mass=bjet_v4_all[3].M()
      j4_b_pt=bjet_v4_all[3].Pt()
      j4_b_eta=bjet_v4_all[3].Eta()
      j4_b_phi=bjet_v4_all[3].Phi()
      mj1j2_b=(bjet_v4_all[0]+bjet_v4_all[1]).M()
    elif (n_bjet_DeepB_M+n_bjet_DeepB_L)==3:
      j1_b_pt=bjet_v4_all[0].Pt()
      j1_b_eta=bjet_v4_all[0].Eta()
      j1_b_phi=bjet_v4_all[0].Phi()
      j1_b_mass=bjet_v4_all[0].M()
      j2_b_pt=bjet_v4_all[1].Pt()
      j2_b_eta=bjet_v4_all[1].Eta()
      j2_b_phi=bjet_v4_all[1].Phi()
      j3_b_mass=bjet_v4_all[2].M()
      j3_b_pt=bjet_v4_all[2].Pt()
      j3_b_eta=bjet_v4_all[2].Eta()
      j3_b_phi=bjet_v4_all[2].Phi()
      mj1j2_b=(bjet_v4_all[0]+bjet_v4_all[1]).M()
    elif (n_bjet_DeepB_M+n_bjet_DeepB_L)==2:
      j1_b_pt=bjet_v4_all[0].Pt()
      j1_b_eta=bjet_v4_all[0].Eta()
      j1_b_phi=bjet_v4_all[0].Phi()
      j1_b_mass=bjet_v4_all[0].M()
      j2_b_pt=bjet_v4_all[1].Pt()
      j2_b_eta=bjet_v4_all[1].Eta()
      j2_b_phi=bjet_v4_all[1].Phi()
      mj1j2_b=(bjet_v4_all[0]+bjet_v4_all[1]).M()
    elif (n_bjet_DeepB_M+n_bjet_DeepB_L)==1:
      j1_b_pt=bjet_v4_all[0].Pt()
      j1_b_eta=bjet_v4_all[0].Eta()
      j1_b_phi=bjet_v4_all[0].Phi()
      j1_b_mass=bjet_v4_all[0].M()
      
      
    if n_tight_nob>3:
      j1_nob_pt=event.Jet_pt_nom[tightJets_nob_id[0]]
      j1_nob_eta=event.Jet_eta[tightJets_nob_id[0]]
      j1_nob_phi=event.Jet_phi[tightJets_nob_id[0]]
      j1_nob_mass=event.Jet_mass_nom[tightJets_nob_id[0]]
      j2_nob_pt=event.Jet_pt_nom[tightJets_nob_id[1]]
      j2_nob_eta=event.Jet_eta[tightJets_nob_id[1]]
      j2_nob_phi=event.Jet_phi[tightJets_nob_id[1]]
      j2_nob_mass=event.Jet_mass_nom[tightJets_nob_id[1]]
      j3_nob_pt=event.Jet_pt_nom[tightJets_nob_id[2]]
      j3_nob_eta=event.Jet_eta[tightJets_nob_id[2]]
      j3_nob_phi=event.Jet_phi[tightJets_nob_id[2]]
      j3_nob_mass=event.Jet_mass_nom[tightJets_nob_id[2]]
      j4_nob_pt=event.Jet_pt_nom[tightJets_nob_id[3]]
      j4_nob_eta=event.Jet_eta[tightJets_nob_id[3]]
      j4_nob_phi=event.Jet_phi[tightJets_nob_id[3]]
      j4_nob_mass=event.Jet_mass_nom[tightJets_nob_id[3]]
      mj1j2_nob=(nobjet_v4_all[0]+nobjet_v4_all[1]).M()
    elif n_tight_nob==3:
      j1_nob_pt=event.Jet_pt_nom[tightJets_nob_id[0]]
      j1_nob_eta=event.Jet_eta[tightJets_nob_id[0]]
      j1_nob_phi=event.Jet_phi[tightJets_nob_id[0]]
      j1_nob_mass=event.Jet_mass_nom[tightJets_nob_id[0]]
      j2_nob_pt=event.Jet_pt_nom[tightJets_nob_id[1]]
      j2_nob_eta=event.Jet_eta[tightJets_nob_id[1]]
      j2_nob_phi=event.Jet_phi[tightJets_nob_id[1]]
      j2_nob_mass=event.Jet_mass_nom[tightJets_nob_id[1]]
      j3_nob_pt=event.Jet_pt_nom[tightJets_nob_id[2]]
      j3_nob_eta=event.Jet_eta[tightJets_nob_id[2]]
      j3_nob_phi=event.Jet_phi[tightJets_nob_id[2]]
      j3_nob_mass=event.Jet_mass_nom[tightJets_nob_id[2]]
      mj1j2_nob=(nobjet_v4_all[0]+nobjet_v4_all[1]).M()
    elif n_tight_nob==2:
      j1_nob_pt=event.Jet_pt_nom[tightJets_nob_id[0]]
      j1_nob_eta=event.Jet_eta[tightJets_nob_id[0]]
      j1_nob_phi=event.Jet_phi[tightJets_nob_id[0]]
      j1_nob_mass=event.Jet_mass_nom[tightJets_nob_id[0]]
      j2_nob_pt=event.Jet_pt_nom[tightJets_nob_id[1]]
      j2_nob_eta=event.Jet_eta[tightJets_nob_id[1]]
      j2_nob_phi=event.Jet_phi[tightJets_nob_id[1]]
      j2_nob_mass=event.Jet_mass_nom[n_tight_nob[1]]
      mj1j2_nob=(nobjet_v4_all[0]+nobjet_v4_all[1]).M()
    elif n_tight_nob==1:
      j1_nob_pt=event.Jet_pt_nom[tightJets_nob_id[0]]
      j1_nob_eta=event.Jet_eta[tightJets_nob_id[0]]
      j1_nob_phi=event.Jet_phi[tightJets_nob_id[0]]
      j1_nob_mass=event.Jet_mass_nom[tightJets_nob_id[0]]

    self.out.fillBranch("j1_nob_pt",j1_nob_pt)
    self.out.fillBranch("j1_nob_eta",j1_nob_eta)
    self.out.fillBranch("j1_nob_phi",j1_nob_phi)
    self.out.fillBranch("j1_nob_mass",j1_nob_mass)
    self.out.fillBranch("j2_nob_pt",j2_nob_pt)
    self.out.fillBranch("j2_nob_eta",j2_nob_eta)
    self.out.fillBranch("j2_nob_phi",j2_nob_phi)
    self.out.fillBranch("j2_nob_mass",j2_nob_mass)
    self.out.fillBranch("j3_nob_pt",j3_nob_pt)
    self.out.fillBranch("j3_nob_eta",j3_nob_eta)
    self.out.fillBranch("j3_nob_phi",j3_nob_phi)
    self.out.fillBranch("j3_nob_mass",j3_nob_mass)
    self.out.fillBranch("j4_nob_pt",j4_nob_pt)
    self.out.fillBranch("j4_nob_eta",j4_nob_eta)
    self.out.fillBranch("j4_nob_phi",j4_nob_phi)
    self.out.fillBranch("j4_nob_mass",j4_nob_mass)
    self.out.fillBranch("mj1j2_nob",mj1j2_nob)

    self.out.fillBranch("j1_b_pt",j1_b_pt)
    self.out.fillBranch("j1_b_eta",j1_b_eta)
    self.out.fillBranch("j1_b_phi",j1_b_phi)
    self.out.fillBranch("j1_b_mass",j1_b_mass)
    self.out.fillBranch("j2_b_pt",j2_b_pt)
    self.out.fillBranch("j2_b_eta",j2_b_eta)
    self.out.fillBranch("j2_b_phi",j2_b_phi)
    self.out.fillBranch("j2_b_mass",j2_b_mass)
    self.out.fillBranch("j3_b_pt",j3_b_pt)
    self.out.fillBranch("j3_b_eta",j3_b_eta)
    self.out.fillBranch("j3_b_phi",j3_b_phi)
    self.out.fillBranch("j3_b_mass",j3_b_mass)
    self.out.fillBranch("j4_b_pt",j4_b_pt)
    self.out.fillBranch("j4_b_eta",j4_b_eta)
    self.out.fillBranch("j4_b_phi",j4_b_phi)
    self.out.fillBranch("j4_b_mass",j4_b_mass)
    self.out.fillBranch("mj1j2_b",mj1j2_b)


    tightJets_id.extend(np.zeros(event.nJet-len(tightJets_id),int)-1)
    tightJets_nob_id.extend(np.zeros(event.nJet-len(tightJets_nob_id),int)-1)
    tightJets_b_DeepCSVmedium_id.extend(np.zeros(event.nJet-len(tightJets_b_DeepCSVmedium_id),int)-1)
    tightJets_b_DeepCSVloose_id.extend(np.zeros(event.nJet-len(tightJets_b_DeepCSVloose_id),int)-1)

    self.out.fillBranch("tightJets_id",tightJets_id)
    self.out.fillBranch("tightJets_nob_id",tightJets_nob_id)
    self.out.fillBranch("tightJets_b_DeepCSVmedium_id",tightJets_b_DeepCSVmedium_id)
    self.out.fillBranch("tightJets_b_DeepCSVloose_id",tightJets_b_DeepCSVloose_id)

    if len(tightLeptons)+len(fakeableLeptons)<1:return False

    Onelep_region=False
    Twolep_region=False
    fourJet_region=False

    Onelep_l1_faketag=False
    Onelep_l1_pt=-99
    Onelep_l1_eta=-99
    Onelep_l1_phi=-99
    Onelep_l1_mass=-99
    Onelep_l1_id=-1
    Onelep_l1_pdgid=-99
    Onelep_l1_isprompt=0

    # Twolep_channel, "1" means two muons, "2" means 1 muon + 1 electron, "3" means two electrons
    Twolep_channel=0
    Twolep_2P0F=False
    Twolep_1P1F=False
    Twolep_0P2F=False
    Twolep_l1_faketag=False
    Twolep_l1_pt=-99
    Twolep_l1_eta=-99
    Twolep_l1_phi=-99
    Twolep_l1_mass=-99
    Twolep_l1_id=-1
    Twolep_l1_pdgid=-99
    Twolep_l1_isprompt=0
    Twolep_l2_faketag=False
    Twolep_l2_pt=-99
    Twolep_l2_eta=-99
    Twolep_l2_phi=-99
    Twolep_l2_mass=-99
    Twolep_l2_id=-1
    Twolep_l2_pdgid=-99
    Twolep_l2_isprompt=0
    Twolep_mll=-99

    # only one tight lepton with pt>30, and no other lepton
    if (sum( lep.Pt() > 30 for lep in tightLeptons )+sum( lep.Pt() > 30 for lep in fakeableLeptons ))==1 and sum( lep.Pt() < 30 for lep in tightLeptons )==0 and sum( lep.Pt() < 30 for lep in fakeableLeptons )==0 and len(looseLeptons)==0:
      Onelep_region=True
      if sum( lep.Pt() > 30 for lep in tightLeptons )==1:
        Onelep_l1_faketag=True
        Onelep_l1_pt=tightLeptons[0].Pt()
        Onelep_l1_eta=tightLeptons[0].Eta()
        Onelep_l1_phi=tightLeptons[0].Phi()
        Onelep_l1_mass=tightLeptons[0].M()
        if sum( lep.Pt() > 30 for lep in tightMuons )==1:
          Onelep_l1_id=tightMuons_id[0]
          Onelep_l1_pdgid=tightMuons_pdgid[0]
          if self.is_mc and event.Muon_genPartFlav[Onelep_l1_id]==1 or event.Muon_genPartFlav[Onelep_l1_id]==15:
            Onelep_l1_isprompt=1
        else:
          Onelep_l1_id=tightElectrons_id[0]
          Onelep_l1_pdgid=tightElectrons_pdgid[0]
          if self.is_mc and event.Electron_genPartFlav[Onelep_l1_id]==1 or event.Electron_genPartFlav[Onelep_l1_id]==15:
            Onelep_l1_isprompt=1
      else:
        Onelep_l1_faketag=False
        Onelep_l1_pt=fakeableLeptons[0].Pt()
        Onelep_l1_eta=fakeableLeptons[0].Eta()
        Onelep_l1_phi=fakeableLeptons[0].Phi()
        Onelep_l1_mass=fakeableLeptons[0].M()
        if sum( lep.Pt() > 30 for lep in fakeable_Muons )==1:
          Onelep_l1_id=fakeable_Muons_id[0]
          Onelep_l1_pdgid=fakeable_Muons_pdgid[0]
          if self.is_mc and event.Muon_genPartFlav[Onelep_l1_id]==1 or event.Muon_genPartFlav[Onelep_l1_id]==15:
            Onelep_l1_isprompt=1
        else:
          Onelep_l1_id=fakeable_Electrons_id[0]
          Onelep_l1_pdgid=fakeable_Electrons_pdgid[0]
          if self.is_mc and event.Electron_genPartFlav[Onelep_l1_id]==1 or event.Electron_genPartFlav[Onelep_l1_id]==15:
            Onelep_l1_isprompt=1


    # only two leptons
    if (sum( lep.Pt() > 25 for lep in tightLeptons )+sum( lep.Pt() > 25 for lep in fakeableLeptons ))==2 and sum( lep.Pt() < 25 for lep in tightLeptons )==0 and sum( lep.Pt() < 25 for lep in fakeableLeptons )==0 and len(looseLeptons)==0:
      Twolep_region=True

      if sum( lep.Pt() > 25 for lep in tightLeptons )==2:
        Twolep_2P0F=True
        Twolep_1P1F=False
        Twolep_0P2F=False
        Twolep_l1_faketag=True
        Twolep_l2_faketag=True
        
        # two muons case
        if sum( lep.Pt() > 25 for lep in tightMuons )==2:
          Twolep_channel=1
          Twolep_l1_pt=tightMuons[0].Pt()
          Twolep_l1_eta=tightMuons[0].Eta()
          Twolep_l1_phi=tightMuons[0].Phi()
          Twolep_l1_mass=tightMuons[0].M()
          Twolep_l1_id=tightMuons_id[0]
          Twolep_l1_pdgid=tightMuons_pdgid[0]
          Twolep_l2_pt=tightMuons[1].Pt()
          Twolep_l2_eta=tightMuons[1].Eta()
          Twolep_l2_phi=tightMuons[1].Phi()
          Twolep_l2_mass=tightMuons[1].M()
          Twolep_l2_id=tightMuons_id[1]
          Twolep_l2_pdgid=tightMuons_pdgid[1]
          Twolep_mll=(tightMuons[0]+tightMuons[1]).M()
          if self.is_mc and event.Muon_genPartFlav[Twolep_l1_id]==1 or event.Muon_genPartFlav[Twolep_l1_id]==15:
            Twolep_l1_isprompt=1
          if self.is_mc and event.Muon_genPartFlav[Twolep_l2_id]==1 or event.Muon_genPartFlav[Twolep_l2_id]==15:
            Twolep_l2_isprompt=1

        # one muon case, in this case, muon is always l1
        if sum( lep.Pt() > 25 for lep in tightMuons )==1:
          Twolep_channel=2
          Twolep_l1_pt=tightMuons[0].Pt()
          Twolep_l1_eta=tightMuons[0].Eta()
          Twolep_l1_phi=tightMuons[0].Phi()
          Twolep_l1_mass=tightMuons[0].M()
          Twolep_l1_id=tightMuons_id[0]
          Twolep_l1_pdgid=tightMuons_pdgid[0]
          Twolep_l2_pt=tightElectrons[0].Pt()
          Twolep_l2_eta=tightElectrons[0].Eta()
          Twolep_l2_phi=tightElectrons[0].Phi()
          Twolep_l2_mass=tightElectrons[0].M()
          Twolep_l2_id=tightElectrons_id[0]
          Twolep_l2_pdgid=tightElectrons_pdgid[0]
          Twolep_mll=(tightMuons[0]+tightElectrons[0]).M()
          if self.is_mc and event.Muon_genPartFlav[Twolep_l1_id]==1 or event.Muon_genPartFlav[Twolep_l1_id]==15:
            Twolep_l1_isprompt=1
          if self.is_mc and event.Electron_genPartFlav[Twolep_l2_id]==1 or event.Electron_genPartFlav[Twolep_l2_id]==15:
            Twolep_l2_isprompt=1

        # two electrons case
        if sum( lep.Pt() > 25 for lep in tightElectrons )==2:
          Twolep_channel=3
          Twolep_l1_pt=tightElectrons[0].Pt()
          Twolep_l1_eta=tightElectrons[0].Eta()
          Twolep_l1_phi=tightElectrons[0].Phi()
          Twolep_l1_mass=tightElectrons[0].M()
          Twolep_l1_id=tightElectrons_id[0]
          Twolep_l1_pdgid=tightElectrons_pdgid[0]
          Twolep_l2_pt=tightElectrons[1].Pt()
          Twolep_l2_eta=tightElectrons[1].Eta()
          Twolep_l2_phi=tightElectrons[1].Phi()
          Twolep_l2_mass=tightElectrons[1].M()
          Twolep_l2_id=tightElectrons_id[1]
          Twolep_l2_pdgid=tightElectrons_pdgid[1]
          Twolep_mll=(tightElectrons[0]+tightElectrons[1]).M()
          if self.is_mc and event.Electron_genPartFlav[Twolep_l1_id]==1 or event.Electron_genPartFlav[Twolep_l1_id]==15:
            Twolep_l1_isprompt=1
          if self.is_mc and event.Electron_genPartFlav[Twolep_l2_id]==1 or event.Electron_genPartFlav[Twolep_l2_id]==15:
            Twolep_l2_isprompt=1

      # one real lepton and one fake lepton
      if sum( lep.Pt() > 25 for lep in tightLeptons )==1:
        Twolep_2P0F=False
        Twolep_1P1F=True
        Twolep_0P2F=False

        if sum( lep.Pt() > 25 for lep in tightElectrons )==0:

          if sum( lep.Pt() > 25 for lep in fakeable_Electrons )==0:
            # one real muon and one fake muon
            Twolep_channel=1
            if fakeable_Muons[0].Pt()>tightMuons[0].Pt():
              Twolep_l1_faketag=False
              Twolep_l2_faketag=True
              Twolep_l1_pt=fakeable_Muons[0].Pt()
              Twolep_l1_eta=fakeable_Muons[0].Eta()
              Twolep_l1_phi=fakeable_Muons[0].Phi()
              Twolep_l1_mass=fakeable_Muons[0].M()
              Twolep_l1_pdgid=fakeable_Muons_pdgid[0]
              Twolep_l1_id=fakeable_Muons_id[0]
              Twolep_l2_pt=tightMuons[0].Pt()
              Twolep_l2_eta=tightMuons[0].Eta()
              Twolep_l2_phi=tightMuons[0].Phi()
              Twolep_l2_mass=tightMuons[0].M()
              Twolep_l2_pdgid=tightMuons_pdgid[0]
              Twolep_l2_id=tightMuons_id[0]
              Twolep_mll=(tightMuons[0]+fakeable_Muons[0]).M()
              if self.is_mc and event.Muon_genPartFlav[Twolep_l1_id]==1 or event.Muon_genPartFlav[Twolep_l1_id]==15:
                Twolep_l1_isprompt=1
              if self.is_mc and event.Muon_genPartFlav[Twolep_l2_id]==1 or event.Muon_genPartFlav[Twolep_l2_id]==15:
                Twolep_l2_isprompt=1
            else:
              Twolep_l1_isfake=True
              Twolep_l2_isfake=False
              Twolep_l1_pt=tightMuons[0].Pt()
              Twolep_l1_eta=tightMuons[0].Eta()
              Twolep_l1_phi=tightMuons[0].Phi()
              Twolep_l1_mass=tightMuons[0].M()
              Twolep_l1_pdgid=tightMuons_pdgid[0]
              Twolep_l1_id=tightMuons_id[0]
              Twolep_l2_pt=fakeable_Muons[0].Pt()
              Twolep_l2_eta=fakeable_Muons[0].Eta()
              Twolep_l2_phi=fakeable_Muons[0].Phi()
              Twolep_l2_mass=fakeable_Muons[0].M()
              Twolep_l2_pdgid=fakeable_Muons_pdgid[0]
              Twolep_l2_id=fakeable_Muons_id[0]
              Twolep_mll=(tightMuons[0]+fakeable_Muons[0]).M()
              if self.is_mc and event.Muon_genPartFlav[Twolep_l1_id]==1 or event.Muon_genPartFlav[Twolep_l1_id]==15:
                Twolep_l1_isprompt=1
              if self.is_mc and event.Muon_genPartFlav[Twolep_l2_id]==1 or event.Muon_genPartFlav[Twolep_l2_id]==15:
                Twolep_l2_isprompt=1
          

          if sum( lep.Pt() > 25 for lep in fakeable_Electrons )==1:
            # one real muon and one fake electron
            Twolep_channel=2
            Twolep_l1_faketag=True
            Twolep_l2_faketag=False
            Twolep_l1_pt=tightMuons[0].Pt()
            Twolep_l1_eta=tightMuons[0].Eta()
            Twolep_l1_phi=tightMuons[0].Phi()
            Twolep_l1_mass=tightMuons[0].M()
            Twolep_l1_pdgid=tightMuons_pdgid[0]
            Twolep_l1_id=tightMuons_id[0]
            Twolep_l2_pt=fakeable_Electrons[0].Pt()
            Twolep_l2_eta=fakeable_Electrons[0].Eta()
            Twolep_l2_phi=fakeable_Electrons[0].Phi()
            Twolep_l2_mass=fakeable_Electrons[0].M()
            Twolep_l2_pdgid=fakeable_Electrons_pdgid[0]
            Twolep_l2_id=fakeable_Electrons_id[0]
            Twolep_mll=(tightMuons[0]+fakeable_Electrons[0]).M()
            if self.is_mc and event.Muon_genPartFlav[Twolep_l1_id]==1 or event.Muon_genPartFlav[Twolep_l1_id]==15:
              Twolep_l1_isprompt=1
            if self.is_mc and event.Electron_genPartFlav[Twolep_l2_id]==1 or event.Electron_genPartFlav[Twolep_l2_id]==15:
              Twolep_l2_isprompt=1


        # one real electron
        else:

          if sum( lep.Pt() > 25 for lep in fakeable_Electrons )==1:
            # one real electron and one fake electron
            Twolep_channel=3
            if fakeable_Electrons[0].Pt()>tightElectrons[0].Pt():
              Twolep_l1_faketag=False
              Twolep_l2_faketag=True
              Twolep_l1_pt=fakeable_Electrons[0].Pt()
              Twolep_l1_eta=fakeable_Electrons[0].Eta()
              Twolep_l1_phi=fakeable_Electrons[0].Phi()
              Twolep_l1_mass=fakeable_Electrons[0].M()
              Twolep_l1_pdgid=fakeable_Electrons_pdgid[0]
              Twolep_l1_id=fakeable_Electrons_id[0]
              Twolep_l2_pt=tightElectrons[0].Pt()
              Twolep_l2_eta=tightElectrons[0].Eta()
              Twolep_l2_phi=tightElectrons[0].Phi()
              Twolep_l2_mass=tightElectrons[0].M()
              Twolep_l2_pdgid=tightElectrons_pdgid[0]
              Twolep_l2_id=tightElectrons_id[0]
              Twolep_mll=(tightElectrons[0]+fakeable_Electrons[0]).M()
              if self.is_mc and event.Electron_genPartFlav[Twolep_l1_id]==1 or event.Electron_genPartFlav[Twolep_l1_id]==15:
                Twolep_l1_isprompt=1
              if self.is_mc and event.Electron_genPartFlav[Twolep_l2_id]==1 or event.Electron_genPartFlav[Twolep_l2_id]==15:
                Twolep_l2_isprompt=1
            else:
              Twolep_l1_faketag=True
              Twolep_l2_faketag=False
              Twolep_l1_pt=tightElectrons[0].Pt()
              Twolep_l1_eta=tightElectrons[0].Eta()
              Twolep_l1_phi=tightElectrons[0].Phi()
              Twolep_l1_mass=tightElectrons[0].M()
              Twolep_l1_pdgid=tightElectrons_pdgid[0]
              Twolep_l1_id=tightElectrons_id[0]
              Twolep_l2_pt=fakeable_Electrons[0].Pt()
              Twolep_l2_eta=fakeable_Electrons[0].Eta()
              Twolep_l2_phi=fakeable_Electrons[0].Phi()
              Twolep_l2_mass=fakeable_Electrons[0].M()
              Twolep_l2_pdgid=fakeable_Electrons_pdgid[0]
              Twolep_l2_id=fakeable_Electrons_id[0]
              Twolep_mll=(tightElectrons[0]+fakeable_Electrons[0]).M()
              if self.is_mc and event.Electron_genPartFlav[Twolep_l1_id]==1 or event.Electron_genPartFlav[Twolep_l1_id]==15:
                Twolep_l1_isprompt=1
              if self.is_mc and event.Electron_genPartFlav[Twolep_l2_id]==1 or event.Electron_genPartFlav[Twolep_l2_id]==15:
                Twolep_l2_isprompt=1
          

          if sum( lep.Pt() > 25 for lep in fakeable_Electrons )==0:
            # one real electron and one fake muon
            Twolep_channel=2
            Twolep_l1_faketag=False
            Twolep_l2_faketag=True
            Twolep_l1_pt=tightElectrons[0].Pt()
            Twolep_l1_eta=tightElectrons[0].Eta()
            Twolep_l1_phi=tightElectrons[0].Phi()
            Twolep_l1_mass=tightElectrons[0].M()
            Twolep_l1_pdgid=tightElectrons_pdgid[0]
            Twolep_l1_id=tightElectrons_id[0]
            Twolep_l2_pt=fakeable_Muons[0].Pt()
            Twolep_l2_eta=fakeable_Muons[0].Eta()
            Twolep_l2_phi=fakeable_Muons[0].Phi()
            Twolep_l2_mass=fakeable_Muons[0].M()
            Twolep_l2_pdgid=fakeable_Muons_pdgid[0]
            Twolep_l2_id=fakeable_Muons_id[0]
            Twolep_mll=(tightElectrns[0]+fakeable_Muons[0]).M()
            if self.is_mc and event.Muon_genPartFlav[Twolep_l1_id]==1 or event.Muon_genPartFlav[Twolep_l1_id]==15:
              Twolep_l1_isprompt=1
            if self.is_mc and event.Electron_genPartFlav[Twolep_l2_id]==1 or event.Electron_genPartFlav[Twolep_l2_id]==15:
              Twolep_l2_isprompt=1

      # two fake leptons
      if sum( lep.Pt() > 25 for lep in tightLeptons )==0:
        Twolep_2P0F=False
        Twolep_1P1F=False
        Twolep_0P2F=True
        Twolep_l1_faketag=False
        Twolep_l2_faketag=False

        # two muons case
        if sum( lep.Pt() > 25 for lep in fakeable_Muons )==2:
          Twolep_channel=1
          Twolep_l1_pt=fakeable_Muons[0].Pt()
          Twolep_l1_eta=fakeable_Muons[0].Eta()
          Twolep_l1_phi=fakeable_Muons[0].Phi()
          Twolep_l1_mass=fakeable_Muons[0].M()
          Twolep_l1_id=fakeable_Muons_id[0]
          Twolep_l1_pdgid=fakeable_Muons_pdgid[0]
          Twolep_l2_pt=fakeable_Muons[1].Pt()
          Twolep_l2_eta=fakeable_Muons[1].Eta()
          Twolep_l2_phi=fakeable_Muons[1].Phi()
          Twolep_l2_mass=fakeable_Muons[1].M()
          Twolep_l2_id=fakeable_Muons_id[1]
          Twolep_l2_pdgid=fakeable_Muons_pdgid[1]
          Twolep_mll=(fakeable_Muons[0]+fakeable_Muons[1]).M()
          if self.is_mc and event.Muon_genPartFlav[Twolep_l1_id]==1 or event.Muon_genPartFlav[Twolep_l1_id]==15:
            Twolep_l1_isprompt=1
          if self.is_mc and event.Muon_genPartFlav[Twolep_l2_id]==1 or event.Muon_genPartFlav[Twolep_l2_id]==15:
            Twolep_l2_isprompt=1

        # one muon case, in this case, muon is always l1
        if sum( lep.Pt() > 25 for lep in fakeable_Muons )==1:
          Twolep_channel=2
          Twolep_l1_pt=fakeable_Muons[0].Pt()
          Twolep_l1_eta=fakeable_Muons[0].Eta()
          Twolep_l1_phi=fakeable_Muons[0].Phi()
          Twolep_l1_mass=fakeable_Muons[0].M()
          Twolep_l1_id=fakeable_Muons_id[0]
          Twolep_l1_pdgid=fakeable_Muons_pdgid[0]
          Twolep_l2_pt=fakeable_Electrons[0].Pt()
          Twolep_l2_eta=fakeable_Electrons[0].Eta()
          Twolep_l2_phi=fakeable_Electrons[0].Phi()
          Twolep_l2_mass=fakeable_Electrons[0].M()
          Twolep_l2_id=fakeable_Electrons_id[0]
          Twolep_l2_pdgid=fakeable_Electrons_pdgid[0]
          Twolep_mll=(fakeable_Muons[0]+fakeable_Electrons[0]).M()
          if self.is_mc and event.Muon_genPartFlav[Twolep_l1_id]==1 or event.Muon_genPartFlav[Twolep_l1_id]==15:
            Twolep_l1_isprompt=1
          if self.is_mc and event.Electron_genPartFlav[Twolep_l2_id]==1 or event.Electron_genPartFlav[Twolep_l2_id]==15:
            Twolep_l2_isprompt=1

        # two electrons case
        if sum( lep.Pt() > 25 for lep in fakeable_Electrons )==2:
          Twolep_channel=3
          Twolep_l1_pt=fakeable_Electrons[0].Pt()
          Twolep_l1_eta=fakeable_Electrons[0].Eta()
          Twolep_l1_phi=fakeable_Electrons[0].Phi()
          Twolep_l1_mass=fakeable_Electrons[0].M()
          Twolep_l1_id=fakeable_Electrons_id[0]
          Twolep_l1_pdgid=fakeable_Electrons_pdgid[0]
          Twolep_l2_pt=fakeable_Electrons[1].Pt()
          Twolep_l2_eta=fakeable_Electrons[1].Eta()
          Twolep_l2_phi=fakeable_Electrons[1].Phi()
          Twolep_l2_mass=fakeable_Electrons[1].M()
          Twolep_l2_id=fakeable_Electrons_id[1]
          Twolep_l2_pdgid=fakeable_Electrons_pdgid[1]
          Twolep_mll=(fakeable_Electrons[0]+fakeable_Electrons[1]).M()
          if self.is_mc and event.Electron_genPartFlav[Twolep_l1_id]==1 or event.Electron_genPartFlav[Twolep_l1_id]==15:
            Twolep_l1_isprompt=1
          if self.is_mc and event.Electron_genPartFlav[Twolep_l2_id]==1 or event.Electron_genPartFlav[Twolep_l2_id]==15:
            Twolep_l2_isprompt=1

    # at least four jets
    if n_tight_jet>3:
      fourJet_region=True


    self.out.fillBranch("Onelep_region",Onelep_region)
    self.out.fillBranch("Twolep_region",Twolep_region)
    self.out.fillBranch("fourJet_region",fourJet_region)
    self.out.fillBranch("Onelep_l1_faketag",Onelep_l1_faketag)
    self.out.fillBranch("Onelep_l1_pt",Onelep_l1_pt)
    self.out.fillBranch("Onelep_l1_eta",Onelep_l1_eta)
    self.out.fillBranch("Onelep_l1_phi",Onelep_l1_phi)
    self.out.fillBranch("Onelep_l1_mass",Onelep_l1_mass)
    self.out.fillBranch("Onelep_l1_id",Onelep_l1_id)
    self.out.fillBranch("Onelep_l1_pdgid",Onelep_l1_pdgid)
    self.out.fillBranch("Onelep_l1_isprompt",Onelep_l1_isprompt)
    self.out.fillBranch("Twolep_channel",Twolep_channel)
    self.out.fillBranch("Twolep_2P0F",Twolep_2P0F)
    self.out.fillBranch("Twolep_1P1F",Twolep_1P1F)
    self.out.fillBranch("Twolep_0P2F",Twolep_0P2F)
    self.out.fillBranch("Twolep_l1_faketag",Twolep_l1_faketag)
    self.out.fillBranch("Twolep_l1_pt",Twolep_l1_pt)
    self.out.fillBranch("Twolep_l1_eta",Twolep_l1_eta)
    self.out.fillBranch("Twolep_l1_phi",Twolep_l1_phi)
    self.out.fillBranch("Twolep_l1_mass",Twolep_l1_mass)
    self.out.fillBranch("Twolep_l1_id",Twolep_l1_id)
    self.out.fillBranch("Twolep_l1_pdgid",Twolep_l1_pdgid)
    self.out.fillBranch("Twolep_l1_isprompt",Twolep_l1_isprompt)
    self.out.fillBranch("Twolep_l2_faketag",Twolep_l2_faketag)
    self.out.fillBranch("Twolep_l2_pt",Twolep_l2_pt)
    self.out.fillBranch("Twolep_l2_eta",Twolep_l2_eta)
    self.out.fillBranch("Twolep_l2_phi",Twolep_l2_phi)
    self.out.fillBranch("Twolep_l2_mass",Twolep_l2_mass)
    self.out.fillBranch("Twolep_l2_id",Twolep_l2_id)
    self.out.fillBranch("Twolep_l2_pdgid",Twolep_l2_pdgid)
    self.out.fillBranch("Twolep_l2_isprompt",Twolep_l2_isprompt)
    self.out.fillBranch("Twolep_mll",Twolep_mll)


    tb_region=False
    tt1L_region=False
    wjet_region=False
    wv_region=False
    tt2L_region=False
    dy_region=False
    met_user=-99
    met_phi_user=-99

    # region definition
    if Onelep_region and fourJet_region and (n_bjet_DeepB_M>0 or n_bjet_DeepB_L>1):
      if n_tight_nob>1:
        if mj1j2_nob>150:
          # signal region, VBS jets invariant mass is large
          tb_region=True
        if mj1j2_nob>60 and mj1j2_nob<120:
          # TT to 1L region, with one of the W boson decays hadronically
          tt1L_region=True

    if Onelep_region and fourJet_region and n_bjet_DeepB_M==0 and n_bjet_DeepB_L==0:
      if n_tight_nob>1:
        # wjets region, no b jets
        wjet_region=True
        if mj1j2_nob>60 and mj1j2_nob<120:
          # WZ/W region, with Z boson decays to two none-b jets
          wv_region=True


    if Twolep_region and fourJet_region and (n_bjet_DeepB_M>0 or n_bjet_DeepB_L>1):
      tt2L_region=True

    if Twolep_region and fourJet_region and n_bjet_DeepB_M==0 and n_bjet_DeepB_L==0 and Twolep_mll>70 and Twolep_mll<110:
      dy_region=True
      

    if self.is_mc:
      met_user=event.MET_T1Smear_pt
      met_phi_user=event.MET_T1Smear_phi
    else:
      met_user=event.MET_T1_pt
      met_phi_user=event.MET_T1_phi

    self.out.fillBranch("tb_region",tb_region)
    self.out.fillBranch("tt1L_region",tt1L_region)
    self.out.fillBranch("wjet_region",wjet_region)
    self.out.fillBranch("wv_region",wv_region)
    self.out.fillBranch("tt2L_region",tt2L_region)
    self.out.fillBranch("dy_region",dy_region)
    self.out.fillBranch("met_user",met_user)
    self.out.fillBranch("met_phi_user",met_phi_user)

    if not (tb_region or tt1L_region or wjet_region or wv_region or tt2L_region or dy_region):
      return False

    return True

TB2016apv = lambda: TBProducer("2016apv")
TB2016 = lambda: TBProducer("2016")
TB2017 = lambda: TBProducer("2017")
TB2018 = lambda: TBProducer("2018")
