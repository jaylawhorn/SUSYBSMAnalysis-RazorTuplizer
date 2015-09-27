// -*- C++ -*-
// Class:      RazorTuplizer
/*
Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Caltech razor team
//         Created:  Thu, 17 Jul 2014 15:00:06 GMT

#ifndef RAZORTUPLIZER_H
#define RAZORTUPLIZER_H

// system include files
#include <memory>
#include <string>
#include <vector>
#include <tuple>

using namespace std;

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//CMSSW package includes
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
//#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
//#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
//#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "SUSYBSMAnalysis/RazorTuplizer/interface/EGammaMvaEleEstimatorCSA14.h"
#include "SUSYBSMAnalysis/RazorTuplizer/interface/ElectronMVAEstimatorRun2NonTrig.h"
//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "SUSYBSMAnalysis/RazorTuplizer/interface/EGammaMvaPhotonEstimator.h"

//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

//------ Array Size Constants ------//
#define OBJECTARRAYSIZE 99
#define GENPARTICLEARRAYSIZE 500

class RazorTuplizer : public edm::EDAnalyzer {
 public:
  //analyzer constructor and destructor
  explicit RazorTuplizer(const edm::ParameterSet&);
  ~RazorTuplizer();

  void loadEvent(const edm::Event& iEvent); //call at the beginning of each event to get input handles from the python config
  virtual void resetBranches();

  //enable desired output variables
  virtual void setBranches();
  virtual void enableEventInfoBranches();
  virtual void enablePileUpBranches();
  virtual void enableMuonBranches();
  virtual void enableElectronBranches();
  virtual void enableTauBranches();
  virtual void enableIsoPFCandidateBranches();
  virtual void enablePhotonBranches();
  virtual void enableJetBranches();
  virtual void enableJetAK8Branches();
  virtual void enableMetBranches();
  virtual void enableTriggerBranches();
  virtual void enableMCBranches();
  virtual void enableGenParticleBranches();

  //select objects and fill tree branches
  virtual bool fillEventInfo(const edm::Event& iEvent);
  virtual bool fillPileUp();//Fill summary PU info
  virtual bool fillMuons();//Fills looseID muon 4-momentum only. PT > 5GeV
  virtual bool fillElectrons();//Fills Ele 4-momentum only. PT > 5GeV
  virtual bool fillTaus();//Fills Tau 4-momentum only. PT > 20GeV
  virtual bool fillIsoPFCandidates();//Fills Isolated PF Candidates, PT > 5 GeV
  virtual bool fillPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup);//Fills photon 4-momentum only. PT > 20GeV && ISO < 0.3
  virtual bool fillJets();//Fills AK5 Jet 4-momentum, CSV, and CISV. PT > 20GeV
  virtual bool fillJetsAK8();//Fills AK8 Jet 4-momentum
  virtual bool fillMet(const edm::Event& iEvent);//Fills MET(mag, phi)
  virtual bool fillTrigger(const edm::Event& iEvent);//Fills trigger information
  virtual bool fillMC();
  virtual bool fillGenParticles();

 protected:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //MVAs for triggering and non-triggering electron ID
  EGammaMvaEleEstimatorCSA14* myMVATrig;
  ElectronMVAEstimatorRun2NonTrig* myMVANonTrig;
  EGammaMvaPhotonEstimator* myPhotonMVA;

  // Control Switches
  bool    isData_;
  bool    useGen_;
  bool enableTriggerInfo_;

  // Mapping of the HLT Triggers and Filters
  string triggerPathNamesFile_;
  string eleHLTFilterNamesFile_;
  string muonHLTFilterNamesFile_;
  string photonHLTFilterNamesFile_;
  static const int NTriggersMAX = 200;
  string triggerPathNames[NTriggersMAX];
  static const int MAX_ElectronHLTFilters = 100;
  string eleHLTFilterNames[MAX_ElectronHLTFilters];
  static const int MAX_MuonHLTFilters = 100;
  string muonHLTFilterNames[MAX_MuonHLTFilters];
  static const int MAX_PhotonHLTFilters = 100;
  string photonHLTFilterNames[MAX_PhotonHLTFilters];

  //EDM tokens for each miniAOD input object
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  edm::EDGetTokenT<reco::GsfElectronCollection> electronsToken_;
  edm::EDGetTokenT<reco::PFTauCollection> tausToken_;
  edm::EDGetTokenT<reco::PhotonCollection> photonsToken_;
  edm::EDGetTokenT<reco::PFJetCollection> jetsToken_;
  edm::EDGetTokenT<reco::PFJetCollection> jetsPuppiToken_;
  edm::EDGetTokenT<reco::PFJetCollection> jetsAK8Token_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken_;
  edm::EDGetTokenT<reco::PFMETCollection> metToken_;
  edm::EDGetTokenT<reco::PFMETCollection> metNoHFToken_;
  edm::EDGetTokenT<reco::PFMETCollection> metPuppiToken_;
  edm::EDGetTokenT<LHEEventProduct> lheInfoToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;
  edm::EDGetTokenT<double> rhoAllToken_;
  edm::EDGetTokenT<double> rhoFastjetAllToken_;
  edm::EDGetTokenT<double> rhoFastjetAllCaloToken_;
  edm::EDGetTokenT<double> rhoFastjetCentralCaloToken_;
  edm::EDGetTokenT<double> rhoFastjetCentralChargedPileUpToken_;
  edm::EDGetTokenT<double> rhoFastjetCentralNeutralToken_;

  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<reco::MuonCollection> muons;
  edm::Handle<reco::GsfElectronCollection> electrons;
  edm::Handle<reco::PFTauCollection> taus;
  edm::Handle<reco::PhotonCollection> photons;
  edm::Handle<reco::PFJetCollection> jets;
  edm::Handle<reco::PFJetCollection> jetsPuppi;
  edm::Handle<reco::PFJetCollection> jetsAK8;
  edm::Handle<reco::PFCandidateCollection> pfCands;
  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<reco::GenJetCollection> genJets;
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<trigger::TriggerEvent> triggerEvent;
  edm::Handle<reco::PFMETCollection> mets;
  edm::Handle<reco::PFMETCollection> metsNoHF;
  edm::Handle<reco::PFMETCollection> metsPuppi;
  edm::Handle<LHEEventProduct> lheInfo;
  edm::Handle<GenEventInfoProduct> genInfo;
  edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
  edm::Handle<double> rhoAll;
  edm::Handle<double> rhoFastjetAll;
  edm::Handle<double> rhoFastjetAllCalo;
  edm::Handle<double> rhoFastjetCentralCalo;
  edm::Handle<double> rhoFastjetCentralChargedPileUp;
  edm::Handle<double> rhoFastjetCentralNeutral;
  const reco::Vertex *myPV;

  //output tree
  TTree *RazorEvents;
  TH1F *NEvents;

  //event info
  bool isData;
  int nPV;
  uint runNum;
  uint lumiNum;
  uint eventNum;
  float pvX;
  float pvY;
  float pvZ;
  float fixedGridRhoAll;
  float fixedGridRhoFastjetAll;
  float fixedGridRhoFastjetAllCalo;
  float fixedGridRhoFastjetCentralCalo;
  float fixedGridRhoFastjetCentralChargedPileUp;
  float fixedGridRhoFastjetCentralNeutral;

  //PU
  int nBunchXing;
  int BunchXing[OBJECTARRAYSIZE];
  int nPU[OBJECTARRAYSIZE];
  float nPUmean[OBJECTARRAYSIZE];

  //Muons
  int nMuons;
  float muonE[OBJECTARRAYSIZE];
  float muonPt[OBJECTARRAYSIZE];
  float muonEta[OBJECTARRAYSIZE];
  float muonPhi[OBJECTARRAYSIZE];
  int muonCharge[OBJECTARRAYSIZE];//muon charge
  bool muonIsLoose[OBJECTARRAYSIZE];
  bool muonIsMedium[OBJECTARRAYSIZE];
  bool muonIsTight[OBJECTARRAYSIZE];
  float muon_d0[OBJECTARRAYSIZE];//transverse impact paramenter
  float muon_dZ[OBJECTARRAYSIZE];//impact parameter
  float muon_ip3d[OBJECTARRAYSIZE];//3d impact paramenter
  float muon_ip3dSignificance[OBJECTARRAYSIZE];//3d impact paramenter/error
  unsigned int muonType[OBJECTARRAYSIZE];//muonTypeBit: global, tracker, standalone
  unsigned int muonQuality[OBJECTARRAYSIZE];//muonID Quality Bits
  float muon_pileupIso[OBJECTARRAYSIZE];
  float muon_chargedIso[OBJECTARRAYSIZE];
  float muon_photonIso[OBJECTARRAYSIZE];
  float muon_neutralHadIso[OBJECTARRAYSIZE];
  float muon_ptrel[OBJECTARRAYSIZE];
  float muon_chargedMiniIso[OBJECTARRAYSIZE];
  float muon_photonAndNeutralHadronMiniIso[OBJECTARRAYSIZE];
  float muon_chargedPileupMiniIso[OBJECTARRAYSIZE];
  float muon_activityMiniIsoAnnulus[OBJECTARRAYSIZE];
  bool  muon_passSingleMuTagFilter[OBJECTARRAYSIZE];
  bool  muon_passHLTFilter[OBJECTARRAYSIZE][MAX_MuonHLTFilters];

  //Electrons
  int nElectrons;
  float eleE[OBJECTARRAYSIZE];
  float elePt[OBJECTARRAYSIZE];
  float eleEta[OBJECTARRAYSIZE];
  float elePhi[OBJECTARRAYSIZE];
  float eleCharge[OBJECTARRAYSIZE];
  float eleE_SC[OBJECTARRAYSIZE];
  //float SC_ElePt[OBJECTARRAYSIZE];
  float eleEta_SC[OBJECTARRAYSIZE];
  float elePhi_SC[OBJECTARRAYSIZE];
  float eleSigmaIetaIeta[OBJECTARRAYSIZE];
  float eleFull5x5SigmaIetaIeta[OBJECTARRAYSIZE];
  float eleR9[OBJECTARRAYSIZE];
  float ele_dEta[OBJECTARRAYSIZE];
  float ele_dPhi[OBJECTARRAYSIZE];
  float ele_HoverE[OBJECTARRAYSIZE];
  float ele_d0[OBJECTARRAYSIZE];
  float ele_dZ[OBJECTARRAYSIZE];
  float ele_ip3d[OBJECTARRAYSIZE];
  float ele_ip3dSignificance[OBJECTARRAYSIZE];
  float ele_pileupIso[OBJECTARRAYSIZE];
  float ele_chargedIso[OBJECTARRAYSIZE];
  float ele_photonIso[OBJECTARRAYSIZE];
  float ele_neutralHadIso[OBJECTARRAYSIZE];
  int ele_MissHits[OBJECTARRAYSIZE];
  bool ele_PassConvVeto[OBJECTARRAYSIZE];
  float ele_OneOverEminusOneOverP[OBJECTARRAYSIZE];
  float ele_IDMVATrig[OBJECTARRAYSIZE];
  float ele_IDMVANonTrig[OBJECTARRAYSIZE];
  float ele_RegressionE[OBJECTARRAYSIZE];
  float ele_CombineP4[OBJECTARRAYSIZE];
  float ele_ptrel[OBJECTARRAYSIZE];
  float ele_chargedMiniIso[OBJECTARRAYSIZE];
  float ele_photonAndNeutralHadronMiniIso[OBJECTARRAYSIZE];
  float ele_chargedPileupMiniIso[OBJECTARRAYSIZE];
  float ele_activityMiniIsoAnnulus[OBJECTARRAYSIZE];
  bool ele_passSingleEleTagFilter[OBJECTARRAYSIZE];
  bool ele_passTPOneTagFilter[OBJECTARRAYSIZE];
  bool ele_passTPTwoTagFilter[OBJECTARRAYSIZE];
  bool ele_passTPOneProbeFilter[OBJECTARRAYSIZE];
  bool ele_passTPTwoProbeFilter[OBJECTARRAYSIZE];
  bool ele_passHLTFilter[OBJECTARRAYSIZE][MAX_ElectronHLTFilters];

  //Taus
  int nTaus;
  float tauE[OBJECTARRAYSIZE];
  float tauPt[OBJECTARRAYSIZE];
  float tauEta[OBJECTARRAYSIZE];
  float tauPhi[OBJECTARRAYSIZE];
  bool tau_IsLoose[OBJECTARRAYSIZE];
  bool tau_IsMedium[OBJECTARRAYSIZE];
  bool tau_IsTight[OBJECTARRAYSIZE];
  bool tau_passEleVetoLoose[OBJECTARRAYSIZE];
  bool tau_passEleVetoMedium[OBJECTARRAYSIZE];
  bool tau_passEleVetoTight[OBJECTARRAYSIZE];
  bool tau_passMuVetoLoose[OBJECTARRAYSIZE];
  bool tau_passMuVetoMedium[OBJECTARRAYSIZE];
  bool tau_passMuVetoTight[OBJECTARRAYSIZE];
  UInt_t tau_ID[OBJECTARRAYSIZE];//tauID Bits
  float tau_combinedIsoDeltaBetaCorr3Hits[OBJECTARRAYSIZE];
  float tau_chargedIsoPtSum[OBJECTARRAYSIZE];
  float tau_neutralIsoPtSum[OBJECTARRAYSIZE];
  float tau_puCorrPtSum[OBJECTARRAYSIZE];
  float tau_eleVetoMVA[OBJECTARRAYSIZE];
  int tau_eleVetoCategory[OBJECTARRAYSIZE];
  float tau_muonVetoMVA[OBJECTARRAYSIZE];
  float tau_isoMVAnewDMwLT[OBJECTARRAYSIZE];
  float tau_isoMVAnewDMwoLT[OBJECTARRAYSIZE];
  float tau_leadCandPt[OBJECTARRAYSIZE];
  int tau_leadCandID[OBJECTARRAYSIZE];
  float tau_leadChargedHadrCandPt[OBJECTARRAYSIZE];
  int tau_leadChargedHadrCandID[OBJECTARRAYSIZE];

  //IsolatedChargedPFCandidates
  int nIsoPFCandidates;
  float isoPFCandidatePt[OBJECTARRAYSIZE];
  float isoPFCandidateEta[OBJECTARRAYSIZE];
  float isoPFCandidatePhi[OBJECTARRAYSIZE];
  float isoPFCandidateIso04[OBJECTARRAYSIZE];
  float isoPFCandidateD0[OBJECTARRAYSIZE];
  int   isoPFCandidatePdgId[OBJECTARRAYSIZE];

  //Photons
  int nPhotons;
  float phoE[OBJECTARRAYSIZE];
  float phoPt[OBJECTARRAYSIZE];
  float phoEta[OBJECTARRAYSIZE];
  float phoPhi[OBJECTARRAYSIZE];
  float phoSigmaIetaIeta[OBJECTARRAYSIZE];
  float phoFull5x5SigmaIetaIeta[OBJECTARRAYSIZE];
  float phoR9[OBJECTARRAYSIZE];
  float pho_HoverE[OBJECTARRAYSIZE];
  float pho_sumChargedHadronPt[OBJECTARRAYSIZE];
  float pho_sumNeutralHadronEt[OBJECTARRAYSIZE];
  float pho_sumPhotonEt[OBJECTARRAYSIZE];
  float pho_sumWorstVertexChargedHadronPt[OBJECTARRAYSIZE];
  bool  pho_isConversion[OBJECTARRAYSIZE];
  bool  pho_passEleVeto[OBJECTARRAYSIZE];
  float pho_RegressionE[OBJECTARRAYSIZE];
  float pho_RegressionEUncertainty[OBJECTARRAYSIZE];
  float pho_IDMVA[OBJECTARRAYSIZE];
  float pho_superClusterEta[OBJECTARRAYSIZE];
  float pho_superClusterPhi[OBJECTARRAYSIZE];
  bool pho_hasPixelSeed[OBJECTARRAYSIZE];
  bool pho_passHLTFilter[OBJECTARRAYSIZE][MAX_PhotonHLTFilters];

  //AK4 Jets
  int nJets;
  float jetE[OBJECTARRAYSIZE];
  float jetPt[OBJECTARRAYSIZE];
  float jetEta[OBJECTARRAYSIZE];
  float jetPhi[OBJECTARRAYSIZE];
  float jetCSV[OBJECTARRAYSIZE];
  float jetCISV[OBJECTARRAYSIZE];
  float jetMass[OBJECTARRAYSIZE];
  float jetJetArea[OBJECTARRAYSIZE];
  float jetPileupE[OBJECTARRAYSIZE];
  float jetPileupId[OBJECTARRAYSIZE];
  int   jetPileupIdFlag[OBJECTARRAYSIZE];
  bool  jetPassIDLoose[OBJECTARRAYSIZE];
  bool  jetPassIDTight[OBJECTARRAYSIZE];
  bool  jetPassMuFrac[OBJECTARRAYSIZE];
  bool  jetPassEleFrac[OBJECTARRAYSIZE];
  int   jetPartonFlavor[OBJECTARRAYSIZE];
  int   jetHadronFlavor[OBJECTARRAYSIZE];
  float jetChargedEMEnergyFraction[OBJECTARRAYSIZE];
  float jetNeutralEMEnergyFraction[OBJECTARRAYSIZE];
  float jetChargedHadronEnergyFraction[OBJECTARRAYSIZE];
  float jetNeutralHadronEnergyFraction[OBJECTARRAYSIZE];
  float jetMuonEnergyFraction[OBJECTARRAYSIZE];
  float jetHOEnergyFraction[OBJECTARRAYSIZE];
  float jetHFHadronEnergyFraction[OBJECTARRAYSIZE];
  float jetHFEMEnergyFraction[OBJECTARRAYSIZE];

  //AK8 Jets
  int nFatJets;
  float fatJetE[OBJECTARRAYSIZE];
  float fatJetPt[OBJECTARRAYSIZE];
  float fatJetEta[OBJECTARRAYSIZE];
  float fatJetPhi[OBJECTARRAYSIZE];
  float fatJetTrimmedM[OBJECTARRAYSIZE];
  float fatJetPrunedM[OBJECTARRAYSIZE];
  float fatJetFilteredM[OBJECTARRAYSIZE];
  float fatJetTau1[OBJECTARRAYSIZE];
  float fatJetTau2[OBJECTARRAYSIZE];
  float fatJetTau3[OBJECTARRAYSIZE];

  //MET
  float metPt;
  float metPhi;
  float sumMET;
  float UncMETdpx;
  float UncMETdpy;
  float UncMETdSumEt;
  float metType0Pt;
  float metType0Phi;
  float metType1Pt;
  float metType1Phi;
  float metType0Plus1Pt;
  float metType0Plus1Phi;
  float metPtRecomputed;
  float metPhiRecomputed;
  float metNoHFPt;
  float metNoHFPhi;
  float metPuppiPt;
  float metPuppiPhi;

  bool Flag_HBHENoiseFilter;
  bool Flag_CSCTightHaloFilter;
  bool Flag_hcalLaserEventFilter;
  bool Flag_EcalDeadCellTriggerPrimitiveFilter;
  bool Flag_goodVertices;
  bool Flag_trackingFailureFilter;
  bool Flag_eeBadScFilter;
  bool Flag_ecalLaserCorrFilter;
  bool Flag_trkPOGFilters;
  bool Flag_trkPOG_manystripclus53X;
  bool Flag_trkPOG_toomanystripclus53X;
  bool Flag_trkPOG_logErrorTooManyClusters;
  bool Flag_METFilters;

  //MC
  int nGenJets;
  float genJetE[OBJECTARRAYSIZE];
  float genJetPt[OBJECTARRAYSIZE];
  float genJetEta[OBJECTARRAYSIZE];
  float genJetPhi[OBJECTARRAYSIZE];
  float genMetPt;
  float genMetPhi;
  float genVertexX;
  float genVertexY;
  float genVertexZ;
  float genWeight;
  unsigned int genSignalProcessID;
  float genQScale;
  float genAlphaQCD;
  float genAlphaQED;

  //gen info
  int nGenParticle;
  int gParticleMotherId[GENPARTICLEARRAYSIZE];
  int gParticleMotherIndex[GENPARTICLEARRAYSIZE];
  int gParticleId[GENPARTICLEARRAYSIZE];
  int gParticleStatus[GENPARTICLEARRAYSIZE];
  float gParticleE[GENPARTICLEARRAYSIZE];
  float gParticlePt[GENPARTICLEARRAYSIZE];
  float gParticleEta[GENPARTICLEARRAYSIZE];
  float gParticlePhi[GENPARTICLEARRAYSIZE];

  //trigger info
  vector<string>  *nameHLT;
  bool triggerDecision[NTriggersMAX];
  int  triggerHLTPrescale[NTriggersMAX];

};

#endif
