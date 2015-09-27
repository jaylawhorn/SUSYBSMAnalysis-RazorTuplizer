#include "RazorTuplizer.h"
//------ Constructors and destructor ------//                                                                                                       
RazorTuplizer::RazorTuplizer(const edm::ParameterSet& iConfig):
  //get inputs from config file
  isData_(iConfig.getParameter<bool> ("isData")),
  useGen_(iConfig.getParameter<bool> ("useGen")),
  verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetsToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetsPuppiToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsPuppi"))),
  jetsAK8Token_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  lheInfoToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheInfo"))),
  genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"))),
  puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfo")))
  //rhoAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoAll"))),
  //rhoFastjetAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAll"))),
  //rhoFastjetAllCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAllCalo"))),
  //rhoFastjetCentralCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralCalo"))),
  //rhoFastjetCentralChargedPileUpToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralChargedPileUp"))),
  //rhoFastjetCentralNeutralToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralNeutral")))
{
  edm::Service<TFileService> fs;

  //set up output tree                                                                                                                              
  RazorEvents = fs->make<TTree>("RazorEvents", "selected AOD information");
  NEvents = fs->make<TH1F>("NEvents",";;NEvents;",1,-0.5,0.5);

}

RazorTuplizer::~RazorTuplizer(){
}

void RazorTuplizer::setBranches(){
  enableEventInfoBranches();
  enablePileUpBranches();
  enableMuonBranches();
  enableElectronBranches();
  enablePhotonBranches();
  enableJetBranches();
  enableJetAK8Branches();

  enableMCBranches();
  enableGenParticleBranches();

}

void RazorTuplizer::enableEventInfoBranches(){
  RazorEvents->Branch("isData", &isData, "isData/O");
  RazorEvents->Branch("nPV", &nPV, "nPV/I");
  RazorEvents->Branch("runNum", &runNum, "runNum/i");
  RazorEvents->Branch("lumiNum", &lumiNum, "lumiNum/i");
  RazorEvents->Branch("eventNum", &eventNum, "eventNum/i");
  RazorEvents->Branch("pvX", &pvX, "pvX/F");
  RazorEvents->Branch("pvY", &pvY, "pvY/F");
  RazorEvents->Branch("pvZ", &pvZ, "pvZ/F");
  RazorEvents->Branch("fixedGridRhoAll", &fixedGridRhoAll, "fixedGridRhoAll/F");
  RazorEvents->Branch("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, "fixedGridRhoFastjetAll/F");
  RazorEvents->Branch("fixedGridRhoFastjetAllCalo", &fixedGridRhoFastjetAllCalo, "fixedGridRhoFastjetAllCalo/F");
  RazorEvents->Branch("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, "fixedGridRhoFastjetCentralCalo/F");
  RazorEvents->Branch("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, "fixedGridRhoFastjetCentralChargedPileUp/F");
  RazorEvents->Branch("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, "fixedGridRhoFastjetCentralNeutral/F");
}                                                                                                                                                 
                                                                                                                                                  
void RazorTuplizer::enablePileUpBranches(){                                                                                                       
  RazorEvents->Branch("nBunchXing", &nBunchXing, "nBunchXing/I");                                                                                 
  RazorEvents->Branch("BunchXing", BunchXing, "BunchXing[nBunchXing]/I");                                                                         
  RazorEvents->Branch("nPU", nPU, "nPU[nBunchXing]/I");                                                                                           
  RazorEvents->Branch("nPUmean", nPUmean, "nPUmean[nBunchXing]/F");                                                                               
}

void RazorTuplizer::enableMuonBranches(){
  RazorEvents->Branch("nMuons", &nMuons,"nMuons/I");
  RazorEvents->Branch("muonE", muonE,"muonE[nMuons]/F");
  RazorEvents->Branch("muonPt", muonPt,"muonPt[nMuons]/F");
  RazorEvents->Branch("muonEta", muonEta,"muonEta[nMuons]/F");
  RazorEvents->Branch("muonPhi", muonPhi,"muonPhi[nMuons]/F");
  RazorEvents->Branch("muonCharge", muonCharge, "muonCharge[nMuons]/I");
  RazorEvents->Branch("muonIsLoose", muonIsLoose,"muonIsLoose[nMuons]/O");
  RazorEvents->Branch("muonIsMedium", muonIsMedium,"muonIsMedium[nMuons]/O");
  RazorEvents->Branch("muonIsTight", muonIsTight,"muonIsTight[nMuons]/O");
  RazorEvents->Branch("muon_d0", muon_d0, "muon_d0[nMuons]/F");
  RazorEvents->Branch("muon_dZ", muon_dZ, "muon_dZ[nMuons]/F");
  RazorEvents->Branch("muon_ip3d", muon_ip3d, "muon_ip3d[nMuons]/F");
  RazorEvents->Branch("muon_ip3dSignificance", muon_ip3dSignificance, "muon_ip3dSignificance[nMuons]/F");
  RazorEvents->Branch("muonType", muonType, "muonType[nMuons]/i");
  RazorEvents->Branch("muonQuality", muonQuality, "muonQuality[nMuons]/i");
  RazorEvents->Branch("muon_pileupIso", muon_pileupIso, "muon_pileupIso[nMuons]/F");
  RazorEvents->Branch("muon_chargedIso", muon_chargedIso, "muon_chargedIso[nMuons]/F");
  RazorEvents->Branch("muon_photonIso", muon_photonIso, "muon_photonIso[nMuons]/F");
  RazorEvents->Branch("muon_neutralHadIso", muon_neutralHadIso, "muon_neutralHadIso[nMuons]/F");
  RazorEvents->Branch("muon_ptrel", muon_ptrel, "muon_ptrel[nMuons]/F");
  RazorEvents->Branch("muon_chargedMiniIso", muon_chargedMiniIso, "muon_chargedMiniIso[nMuons]/F");
  RazorEvents->Branch("muon_photonAndNeutralHadronMiniIso", muon_photonAndNeutralHadronMiniIso, "muon_photonAndNeutralHadronMiniIso[nMuons]/F");
  RazorEvents->Branch("muon_chargedPileupMiniIso", muon_chargedPileupMiniIso, "muon_chargedPileupMiniIso[nMuons]/F");
  RazorEvents->Branch("muon_activityMiniIsoAnnulus", muon_activityMiniIsoAnnulus, "muon_activityMiniIsoAnnulus[nMuons]/F");
  RazorEvents->Branch("muon_passSingleMuTagFilter", muon_passSingleMuTagFilter, "muon_passSingleMuTagFilter[nMuons]/O");
  //RazorEvents->Branch("muon_passHLTFilter", &muon_passHLTFilter, Form("muon_passHLTFilter[nMuons][%d]/O",MAX_MuonHLTFilters));
}

void RazorTuplizer::enableElectronBranches(){
  RazorEvents->Branch("nElectrons", &nElectrons,"nElectrons/I");
  RazorEvents->Branch("eleE", eleE,"eleE[nElectrons]/F");
  RazorEvents->Branch("elePt", elePt,"elePt[nElectrons]/F");
  RazorEvents->Branch("eleEta", eleEta,"eleEta[nElectrons]/F");
  RazorEvents->Branch("elePhi", elePhi,"elePhi[nElectrons]/F");
  RazorEvents->Branch("eleCharge", eleCharge, "eleCharge[nElectrons]/F");
  //RazorEvents->Branch("EleE_SC", eleE_SC,"eleE_SC[nElectrons]/F");                                                                                    
  RazorEvents->Branch("eleEta_SC", eleEta_SC,"eleEta_SC[nElectrons]/F");
  //RazorEvents->Branch("elePhi_SC", elePhi_SC,"elePhi_SC[nElectrons]/F");                                                                              
  RazorEvents->Branch("eleSigmaIetaIeta", eleSigmaIetaIeta, "eleSigmaIetaIeta[nElectrons]/F");
  RazorEvents->Branch("eleFull5x5SigmaIetaIeta", eleFull5x5SigmaIetaIeta, "eleFull5x5SigmaIetaIeta[nElectrons]/F");
  RazorEvents->Branch("eleR9", eleR9, "eleR9[nElectrons]/F");
  RazorEvents->Branch("ele_dEta", ele_dEta, "ele_dEta[nElectrons]/F");
  RazorEvents->Branch("ele_dPhi", ele_dPhi, "ele_dPhi[nElectrons]/F");
  RazorEvents->Branch("ele_HoverE", ele_HoverE, "ele_HoverE[nElectrons]/F");
  RazorEvents->Branch("ele_d0", ele_d0, "ele_d0[nElectrons]/F");
  RazorEvents->Branch("ele_dZ", ele_dZ, "ele_dZ[nElectrons]/F");
  RazorEvents->Branch("ele_ip3d", ele_ip3d, "ele_ip3d[nElectrons]/F");
  RazorEvents->Branch("ele_ip3dSignificance", ele_ip3dSignificance, "ele_ip3dSignificance[nElectrons]/F");
  RazorEvents->Branch("ele_pileupIso", ele_pileupIso, "ele_pileupIso[nElectrons]/F");
  RazorEvents->Branch("ele_chargedIso", ele_chargedIso, "ele_chargedIso[nElectrons]/F");
  RazorEvents->Branch("ele_photonIso", ele_photonIso, "ele_photonIso[nElectrons]/F");
  RazorEvents->Branch("ele_neutralHadIso", ele_neutralHadIso, "ele_neutralHadIso[nElectrons]/F");
  RazorEvents->Branch("ele_MissHits", ele_MissHits, "ele_MissHits[nElectrons]/I");
  RazorEvents->Branch("ele_PassConvVeto", ele_PassConvVeto, "ele_PassConvVeto[nElectrons]/O");
  RazorEvents->Branch("ele_OneOverEminusOneOverP", ele_OneOverEminusOneOverP, "ele_OneOverEminusOneOverP[nElectrons]/F");
  RazorEvents->Branch("ele_IDMVATrig", ele_IDMVATrig, "ele_IDMVATrig[nElectrons]/F");
  RazorEvents->Branch("ele_IDMVANonTrig", ele_IDMVANonTrig, "ele_IDMVANonTrig[nElectrons]/F");
  RazorEvents->Branch("ele_RegressionE", ele_RegressionE, "ele_RegressionE[nElectrons]/F");
  RazorEvents->Branch("ele_CombineP4", ele_CombineP4, "ele_CombineP4[nElectrons]/F");
  RazorEvents->Branch("ele_ptrel", ele_ptrel, "ele_ptrel[nElectrons]/F");
  RazorEvents->Branch("ele_chargedMiniIso", ele_chargedMiniIso, "ele_chargedMiniIso[nElectrons]/F");
  RazorEvents->Branch("ele_photonAndNeutralHadronMiniIso", ele_photonAndNeutralHadronMiniIso, "ele_photonAndNeutralHadronMiniIso[nElectrons]/F");
  RazorEvents->Branch("ele_chargedPileupMiniIso", ele_chargedPileupMiniIso, "ele_chargedPileupMiniIso[nElectrons]/F");
  RazorEvents->Branch("ele_activityMiniIsoAnnulus", ele_activityMiniIsoAnnulus, "ele_activityMiniIsoAnnulus[nElectrons]/F");
  RazorEvents->Branch("ele_passSingleEleTagFilter", ele_passSingleEleTagFilter, "ele_passSingleEleTagFilter[nElectrons]/O");
  RazorEvents->Branch("ele_passTPOneTagFilter", ele_passTPOneTagFilter, "ele_passTPOneTagFilter[nElectrons]/O");
  RazorEvents->Branch("ele_passTPTwoTagFilter", ele_passTPTwoTagFilter, "ele_passTPTwoTagFilter[nElectrons]/O");
  RazorEvents->Branch("ele_passTPOneProbeFilter", ele_passTPOneProbeFilter, "ele_passTPOneProbeFilter[nElectrons]/O");
  RazorEvents->Branch("ele_passTPTwoProbeFilter", ele_passTPTwoProbeFilter, "ele_passTPTwoProbeFilter[nElectrons]/O");
  //RazorEvents->Branch("ele_passHLTFilter", &ele_passHLTFilter, Form("ele_passHLTFilter[nElectrons][%d]/O",MAX_ElectronHLTFilters));
}

void RazorTuplizer::enablePhotonBranches(){
  RazorEvents->Branch("nPhotons", &nPhotons,"nPhotons/I");
  RazorEvents->Branch("phoE", phoE,"phoE[nPhotons]/F");
  RazorEvents->Branch("phoPt", phoPt,"phoPt[nPhotons]/F");
  RazorEvents->Branch("phoEta", phoEta,"phoEta[nPhotons]/F");
  RazorEvents->Branch("phoPhi", phoPhi,"phoPhi[nPhotons]/F");
  RazorEvents->Branch("phoSigmaIetaIeta", phoSigmaIetaIeta, "phoSigmaIetaIeta[nPhotons]/F");
  RazorEvents->Branch("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta, "phoFull5x5SigmaIetaIeta[nPhotons]/F");
  RazorEvents->Branch("phoR9", phoR9, "phoR9[nPhotons]/F");
  RazorEvents->Branch("pho_HoverE", pho_HoverE, "pho_HoverE[nPhotons]/F");
  RazorEvents->Branch("pho_sumChargedHadronPt", pho_sumChargedHadronPt, "pho_sumChargedHadronPt[nPhotons]/F");
  RazorEvents->Branch("pho_sumNeutralHadronEt", pho_sumNeutralHadronEt, "pho_sumNeutralHadronEt[nPhotons]/F");
  RazorEvents->Branch("pho_sumPhotonEt", pho_sumPhotonEt, "pho_sumPhotonEt[nPhotons]/F");
  RazorEvents->Branch("pho_sumWorstVertexChargedHadronPt", pho_sumWorstVertexChargedHadronPt, "pho_sumWorstVertexChargedHadronPt[nPhotons]/F");
  RazorEvents->Branch("pho_isConversion", pho_isConversion, "pho_isConversion[nPhotons]/O");
  RazorEvents->Branch("pho_passEleVeto", pho_passEleVeto, "pho_passEleVeto[nPhotons]/O");
  RazorEvents->Branch("pho_RegressionE", pho_RegressionE, "pho_RegressionE[nPhotons]/F");
  RazorEvents->Branch("pho_RegressionEUncertainty", pho_RegressionEUncertainty, "pho_RegressionEUncertainty[nPhotons]/F");
  RazorEvents->Branch("pho_IDMVA", pho_IDMVA, "pho_IDMVA[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterEta", pho_superClusterEta, "pho_superClusterEta[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterPhi", pho_superClusterPhi, "pho_superClusterPhi[nPhotons]/F");
  RazorEvents->Branch("pho_hasPixelSeed", pho_hasPixelSeed, "pho_hasPixelSeed[nPhotons]/O");
  //RazorEvents->Branch("pho_passHLTFilter", &pho_passHLTFilter, Form("pho_passHLTFilter[nPhotons][%d]/O",MAX_PhotonHLTFilters));
}

void RazorTuplizer::enableJetBranches(){
  RazorEvents->Branch("nJets", &nJets,"nJets/I");
  RazorEvents->Branch("jetE", jetE,"jetE[nJets]/F");
  RazorEvents->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  RazorEvents->Branch("jetEta", jetEta,"jetEta[nJets]/F");
  RazorEvents->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  RazorEvents->Branch("jetCSV", jetCSV,"jetCSV[nJets]/F");
  RazorEvents->Branch("jetCISV", jetCISV,"jetCISV[nJets]/F");
  RazorEvents->Branch("jetMass", jetMass, "jetMass[nJets]/F");
  RazorEvents->Branch("jetJetArea", jetJetArea, "jetJetArea[nJets]/F");
  RazorEvents->Branch("jetPileupE", jetPileupE, "jetPileupE[nJets]/F");
  RazorEvents->Branch("jetPileupId", jetPileupId, "jetPileupId[nJets]/F");
  RazorEvents->Branch("jetPileupIdFlag", jetPileupIdFlag, "jetPileupIdFlag[nJets]/I");
  RazorEvents->Branch("jetPassIDLoose", jetPassIDLoose, "jetPassIDLoose[nJets]/O");
  RazorEvents->Branch("jetPassIDTight", jetPassIDTight, "jetPassIDTight[nJets]/O");
  RazorEvents->Branch("jetPassMuFrac", jetPassMuFrac, "jetPassMuFrac[nJets]/O");
  RazorEvents->Branch("jetPassEleFrac", jetPassEleFrac, "jetPassEleFrac[nJets]/O");
  RazorEvents->Branch("jetPartonFlavor", jetPartonFlavor, "jetPartonFlavor[nJets]/I");
  RazorEvents->Branch("jetHadronFlavor", jetHadronFlavor, "jetHadronFlavor[nJets]/I");
  RazorEvents->Branch("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction, "jetChargedEMEnergyFraction[nJets]/F");
  RazorEvents->Branch("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction, "jetNeutralEMEnergyFraction[nJets]/F");
  RazorEvents->Branch("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction, "jetChargedHadronEnergyFraction[nJets]/F");
  RazorEvents->Branch("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction, "jetNeutralHadronEnergyFraction[nJets]/F");
  RazorEvents->Branch("jetMuonEnergyFraction", jetMuonEnergyFraction, "jetMuonEnergyFraction[nJets]/F");
  RazorEvents->Branch("jetHOEnergyFraction", jetHOEnergyFraction, "jetHOEnergyFraction[nJets]/F");
  RazorEvents->Branch("jetHFHadronEnergyFraction", jetHFHadronEnergyFraction, "jetHFHadronEnergyFraction[nJets]/F");
  RazorEvents->Branch("jetHFEMEnergyFraction",jetHFEMEnergyFraction, "jetHFEMEnergyFraction[nJets]/F");
}

void RazorTuplizer::enableJetAK8Branches(){
  RazorEvents->Branch("nFatJets", &nFatJets,"nFatJets/i");
  RazorEvents->Branch("fatJetE", fatJetE,"fatJetE[nFatJets]/F");
  RazorEvents->Branch("fatJetPt", fatJetPt,"fatJetPt[nFatJets]/F");
  RazorEvents->Branch("fatJetEta", fatJetEta,"fatJetEta[nFatJets]/F");
  RazorEvents->Branch("fatJetPhi", fatJetPhi,"fatJetPhi[nFatJets]/F");
  RazorEvents->Branch("fatJetPrunedM", fatJetPrunedM,"fatJetPrunedM[nFatJets]/F");
  RazorEvents->Branch("fatJetTrimmedM", fatJetTrimmedM,"fatJetTrimmedM[nFatJets]/F");
  RazorEvents->Branch("fatJetFilteredM", fatJetFilteredM,"fatJetFilteredM[nFatJets]/F");
  RazorEvents->Branch("fatJetTau1", fatJetTau1,"fatJetTau1[nFatJets]/F");
  RazorEvents->Branch("fatJetTau2", fatJetTau2,"fatJetTau2[nFatJets]/F");
  RazorEvents->Branch("fatJetTau3", fatJetTau3,"fatJetTau3[nFatJets]/F");
}

void RazorTuplizer::enableMCBranches(){
  RazorEvents->Branch("nGenJets", &nGenJets, "nGenJets/I");
  RazorEvents->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
  RazorEvents->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
  RazorEvents->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
  RazorEvents->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");
  RazorEvents->Branch("genMetPt", &genMetPt, "genMetPt/F");
  RazorEvents->Branch("genMetPhi", &genMetPhi, "genMetPhi/F");
  RazorEvents->Branch("genVertexZ", &genVertexZ, "genVertexZ/F");
  RazorEvents->Branch("genWeight", &genWeight, "genWeight/F");
  RazorEvents->Branch("genSignalProcessID", &genSignalProcessID, "genSignalProcessID/i");
  RazorEvents->Branch("genQScale", &genQScale, "genQScale/F");
  RazorEvents->Branch("genAlphaQCD", &genAlphaQCD, "genAlphaQCD/F");
  RazorEvents->Branch("genAlphaQED", &genAlphaQED, "genAlphaQED/F");
}

void RazorTuplizer::enableGenParticleBranches(){
  RazorEvents->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
  RazorEvents->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  RazorEvents->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  RazorEvents->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  RazorEvents->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  RazorEvents->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  RazorEvents->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  RazorEvents->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");
  RazorEvents->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
}

void RazorTuplizer::loadEvent(const edm::Event& iEvent){
  iEvent.getByToken(verticesToken_, vertices);
  iEvent.getByToken(muonsToken_, muons);
  iEvent.getByToken(electronsToken_, electrons);
  iEvent.getByToken(photonsToken_, photons);
  iEvent.getByToken(jetsToken_, jets);
  iEvent.getByToken(jetsPuppiToken_, jetsPuppi);
  iEvent.getByToken(jetsAK8Token_, jetsAK8);
  if (useGen_) {
    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genJetsToken_,genJets);
    iEvent.getByToken(lheInfoToken_, lheInfo);
    iEvent.getByToken(genInfoToken_,genInfo);
    iEvent.getByToken(puInfoToken_,puInfo);
  }
}

void RazorTuplizer::resetBranches(){
  //Event
  nPV = -1;
  eventNum = 0;
  lumiNum = 0;
  runNum = 0;
  pvX = -99.0;
  pvY = -99.0;
  pvZ = -99.0;
  fixedGridRhoAll = -99.0;
  fixedGridRhoFastjetAll = -99.0;
  fixedGridRhoFastjetAllCalo = -99.0;
  fixedGridRhoFastjetCentralCalo = -99.0;
  fixedGridRhoFastjetCentralChargedPileUp = -99.0;
  fixedGridRhoFastjetCentralNeutral = -99.0;

  nBunchXing = 0;
  nMuons = 0;
  nElectrons = 0;
  nPhotons = 0;
  nJets = 0;
  nFatJets = 0;
  nGenJets = 0;
  nGenParticle = 0;

  for(int i = 0; i < 99; i++){
    //PU
    BunchXing[i] = -99;
    nPU[i] = -99;
    nPUmean[i] = -99.0;

    //Muon
    muonE[i] = 0.0;
    muonPt[i] = 0.0;
    muonEta[i] = 0.0;
    muonPhi[i] = 0.0;
    muonCharge[i] = -99;
    muonIsLoose[i] = false;
    muonIsMedium[i] = false;
    muonIsTight[i] = false;
    muon_d0[i] = -99.0;
    muon_dZ[i] = -99.0;
    muon_ip3d[i] = -99.0;
    muon_ip3dSignificance[i] = -99.0;
    muonType[i] = 0;
    muonQuality[i] = 0;
    muon_pileupIso[i] = -99.0;
    muon_chargedIso[i] = -99.0;
    muon_photonIso[i] = -99.0;
    muon_neutralHadIso[i] = -99.0;
    muon_ptrel[i] = -99.0;
    muon_chargedMiniIso[i] = -99.0;
    muon_photonAndNeutralHadronMiniIso[i] = -99.0;
    muon_chargedPileupMiniIso[i] = -99.0;
    muon_activityMiniIsoAnnulus[i] = -99.0;
    muon_passSingleMuTagFilter[i] = false;
    //for (int q=0;q<MAX_MuonHLTFilters;q++) muon_passHLTFilter[i][q] = false;  

    //Electron                                                                                                                                      
    eleE[i] = 0.0;
    elePt[i] = 0.0;
    eleEta[i] = 0.0;
    elePhi[i] = 0.0;
    eleE_SC[i] = -99.0;
    eleEta_SC[i] = -99.0;
    elePhi_SC[i] = -99.0;
    eleSigmaIetaIeta[i] = -99.0;
    eleFull5x5SigmaIetaIeta[i] = -99.0;
    eleR9[i] = -99;
    ele_dEta[i] = -99;
    ele_dPhi[i] = -99;
    ele_HoverE[i] = -99;
    ele_d0[i] = -99;
    ele_dZ[i] = -99;
    ele_ip3d[i] = -99;
    ele_ip3dSignificance[i] = -99;
    ele_pileupIso[i] = -99.0;
    ele_chargedIso[i] = -99.0;
    ele_photonIso[i] = -99.0;
    ele_neutralHadIso[i] = -99.0;
    ele_MissHits[i] = -99;
    ele_PassConvVeto[i] = false;
    ele_OneOverEminusOneOverP[i] = -99.0;
    ele_IDMVATrig[i] = -99.0;
    ele_IDMVANonTrig[i] = -99.0;
    ele_RegressionE[i] = -99.0;
    ele_CombineP4[i] = -99.0;
    ele_ptrel[i] = -99.0;
    ele_chargedMiniIso[i] = -99.0;
    ele_photonAndNeutralHadronMiniIso[i] = -99.0;
    ele_chargedPileupMiniIso[i] = -99.0;
    ele_activityMiniIsoAnnulus[i] = -99.0;
    ele_passSingleEleTagFilter[i] = false;
    ele_passTPOneTagFilter[i] = false;
    ele_passTPTwoTagFilter[i] = false;
    ele_passTPOneProbeFilter[i] = false;
    ele_passTPTwoProbeFilter[i] = false;
    //for (int q=0;q<MAX_ElectronHLTFilters;q++) ele_passHLTFilter[i][q] = false;

    //Photon                                                                                                                                        
    phoE[i] = 0.0;
    phoPt[i] = 0.0;
    phoEta[i] = 0.0;
    phoPhi[i] = 0.0;
    phoSigmaIetaIeta[i] = -99.0;
    phoFull5x5SigmaIetaIeta[i] = -99.0;
    phoR9[i] = -99.0;
    pho_HoverE[i] = -99.0;
    pho_sumChargedHadronPt[i] = -99.0;
    pho_sumNeutralHadronEt[i] = -99.0;
    pho_sumPhotonEt[i] = -99.0;
    pho_sumWorstVertexChargedHadronPt[i] = -99.0;
    pho_isConversion[i] = false;
    pho_passEleVeto[i] = false;
    pho_RegressionE[i] = -99.0;
    pho_RegressionEUncertainty[i] = -99.0;
    pho_IDMVA[i] = -99.0;
    pho_superClusterEta[i] = -99.0;
    pho_superClusterPhi[i] = -99.0;
    pho_hasPixelSeed[i] = false;
    //for (int q=0;q<MAX_PhotonHLTFilters;q++) pho_passHLTFilter[i][q] = false;

    //Jet                                                                                                                                           
    jetE[i] = 0.0;
    jetPt[i] = 0.0;
    jetEta[i] = 0.0;
    jetPhi[i] = 0.0;
    jetCSV[i] = 0.0;
    jetCISV[i] = 0.0;
    jetMass[i] =  -99.0;
    jetJetArea[i] = -99.0;
    jetPileupE[i] = -99.0;
    jetPileupId[i] = -99.0;
    jetPileupIdFlag[i] = -1;
    jetPassIDLoose[i] = false;
    jetPassIDTight[i] = false;
    jetPassMuFrac[i] = false;
    jetPassEleFrac[i] = false;
    jetPartonFlavor[i] = 0;
    jetHadronFlavor[i] = 0;
    jetChargedEMEnergyFraction[i] = -99.0;
    jetNeutralEMEnergyFraction[i] = -99.0;
    jetChargedHadronEnergyFraction[i] = -99.0;
    jetNeutralHadronEnergyFraction[i] = -99.0;
    jetMuonEnergyFraction[i] = -99.0;
    jetHOEnergyFraction[i] = -99.0;
    jetHFHadronEnergyFraction[i] = -99.0;
    jetHFEMEnergyFraction[i] = -99.0;

    //AK8 Jet                                                                                                                                       
    fatJetE[i] = 0.0;
    fatJetPt[i] = 0.0;
    fatJetEta[i] = 0.0;
    fatJetPhi[i] = 0.0;
    fatJetPrunedM[i] = 0.0;
    fatJetTrimmedM[i] = 0.0;
    fatJetFilteredM[i] = 0.0;
    fatJetTau1[i] = 0.0;
    fatJetTau2[i] = 0.0;
    fatJetTau3[i] = 0.0;

    genJetE[i] = 0.0;
    genJetPt[i] = 0.0;
    genJetEta[i] = 0.0;
    genJetPhi[i] = 0.0;

  }

  for(int i = 0; i < 500; i++){
    //Gen Particle
    gParticleMotherId[i] = -99999;
    gParticleMotherIndex[i] = -99999;
    gParticleId[i] = -99999;
    gParticleStatus[i] = -99999;
    gParticleE[i] = -99999.0;
    gParticlePt[i] = -99999.0;
    gParticleEta[i] = -99999.0;
    gParticlePhi[i] = -99999.0;

  }

}

bool RazorTuplizer::fillEventInfo(const edm::Event& iEvent){
  //store basic event info                                                                                                                          
  isData = isData_;
  runNum = iEvent.id().run();
  lumiNum = iEvent.luminosityBlock();
  eventNum = iEvent.id().event();

  //select the primary vertex, if any                                                                                                               
  nPV = 0;
  bool foundPV = false;
  for(unsigned int i = 0; i < vertices->size(); i++){
    if(vertices->at(i).isValid() && !vertices->at(i).isFake()){
      if (!foundPV) {
        myPV = &(vertices->at(i));
        foundPV = true;
      }
      nPV++;
    }
  }

  if(nPV == 0)return false;
  if (foundPV) {
    pvX = myPV->x();
    pvY = myPV->y();
    pvZ = myPV->z();
  } else {
    return false;
  }

  //get rho                                                                                                                                         
  //fixedGridRhoAll = *rhoAll;
  //fixedGridRhoFastjetAll = *rhoFastjetAll;
  //fixedGridRhoFastjetAllCalo = *rhoFastjetAllCalo;
  //fixedGridRhoFastjetCentralCalo = *rhoFastjetCentralCalo;
  //fixedGridRhoFastjetCentralChargedPileUp = *rhoFastjetCentralChargedPileUp;
  //fixedGridRhoFastjetCentralNeutral = *rhoFastjetCentralNeutral;

  return true;
}

bool RazorTuplizer::fillPileUp(){
  for(const PileupSummaryInfo &pu : *puInfo){
    //std::cout << pu.getBunchCrossing() << std::endl;
    //if (pu.getBunchCrossing() == -1 || pu.getBunchCrossing() == 0 || pu.getBunchCrossing() == 1) {
    BunchXing[nBunchXing] = pu.getBunchCrossing();
    nPU[nBunchXing] = pu.getPU_NumInteractions();
    nPUmean[nBunchXing] = pu.getTrueNumInteractions();
    nBunchXing++;
    //}
  }
  return true;
}

bool RazorTuplizer::fillMuons(){
  for(const reco::Muon &mu : *muons){
    if(mu.pt() < 5) continue;
    muonE[nMuons] = mu.energy();
    muonPt[nMuons] = mu.pt();
    muonEta[nMuons] = mu.eta();
    muonPhi[nMuons] = mu.phi();
    muonCharge[nMuons] = mu.charge();
    //muonIsLoose[nMuons] = mu.isLooseMuon();                                                                                                       
    //muonIsMedium[nMuons] = mu.isMediumMuon();                                                                                                     
    //muonIsTight[nMuons] = mu.isTightMuon(*myPV);                                                                                                  
    //muon_d0[nMuons] = -mu.muonBestTrack()->dxy(myPV->position());
    //muon_dZ[nMuons] = mu.muonBestTrack()->dz(myPV->position());
    //muon_ip3d[nMuons] = mu.dB(pat::Muon::PV3D);                                                                                                   
    //muon_ip3dSignificance[nMuons] = mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D);                                                               
    //muonType[nMuons] = mu.isMuon() + mu.isGlobalMuon() + mu.isTrackerMuon() + mu.isStandAloneMuon()
    //+ mu.isCaloMuon() + mu.isPFMuon() + mu.isRPCMuon();
    //muonQuality[nMuons] =
    //muon::isGoodMuon(mu,muon::All)
    //+ muon::isGoodMuon(mu,muon::AllGlobalMuons)
    //+ muon::isGoodMuon(mu,muon::AllStandAloneMuons)
    //+ muon::isGoodMuon(mu,muon::AllTrackerMuons)
    //+ muon::isGoodMuon(mu,muon::TrackerMuonArbitrated)
    //+ muon::isGoodMuon(mu,muon::AllArbitrated)
    //+ muon::isGoodMuon(mu,muon::GlobalMuonPromptTight)
    //+ muon::isGoodMuon(mu,muon::TMLastStationLoose)
    //+ muon::isGoodMuon(mu,muon::TMLastStationTight)
    //+ muon::isGoodMuon(mu,muon::TM2DCompatibilityLoose)
    //+ muon::isGoodMuon(mu,muon::TM2DCompatibilityTight)
    //+ muon::isGoodMuon(mu,muon::TMOneStationLoose)
    //+ muon::isGoodMuon(mu,muon::TMOneStationTight)
    //+ muon::isGoodMuon(mu,muon::TMLastStationOptimizedLowPtLoose)
    //+ muon::isGoodMuon(mu,muon::TMLastStationOptimizedLowPtTight)
    //+ muon::isGoodMuon(mu,muon::GMTkChiCompatibility)
    //+ muon::isGoodMuon(mu,muon::GMStaChiCompatibility)
    //+ muon::isGoodMuon(mu,muon::GMTkKinkTight)
    //+ muon::isGoodMuon(mu,muon::TMLastStationAngLoose)
    //+ muon::isGoodMuon(mu,muon::TMLastStationAngTight)
    //+ muon::isGoodMuon(mu,muon::TMOneStationAngLoose)
    //+ muon::isGoodMuon(mu,muon::TMOneStationAngTight)
    //+ muon::isGoodMuon(mu,muon::TMLastStationOptimizedBarrelLowPtLoose)
    //+ muon::isGoodMuon(mu,muon::TMLastStationOptimizedBarrelLowPtTight)
    //+ muon::isGoodMuon(mu,muon::RPCMuLoose);
    //muon_pileupIso[nMuons] = mu.pfIsolationR04().sumPUPt;
    //muon_chargedIso[nMuons] = mu.pfIsolationR04().sumChargedHadronPt;
    //muon_photonIso[nMuons] = mu.pfIsolationR04().sumPhotonEt;
    //muon_neutralHadIso[nMuons] = mu.pfIsolationR04().sumNeutralHadronEt;
    //muon_ptrel[nMuons] = getLeptonPtRel( jets, &mu );
    //tuple<double,double,double> PFMiniIso = getPFMiniIsolation(packedPFCands, dynamic_cast<const reco::Candidate *>(&mu), 0.05, 0.2, 10., false, false);
    //muon_chargedMiniIso[nMuons] = std::get<0>(PFMiniIso);
    //muon_photonAndNeutralHadronMiniIso[nMuons] = std::get<1>(PFMiniIso);
    //muon_chargedPileupMiniIso[nMuons] = std::get<2>(PFMiniIso);
    //muon_activityMiniIsoAnnulus[nMuons] = ActivityPFMiniIsolationAnnulus( packedPFCands, dynamic_cast<const reco::Candidate *>(&mu), 0.4, 0.05, 0.2 , 10.);
    nMuons++;
  }
  return true;
}

bool RazorTuplizer::fillElectrons(){
  for(const reco::GsfElectron &ele : *electrons){
    if(ele.pt() < 5) continue;
    eleE[nElectrons] = ele.energy();
    elePt[nElectrons] = ele.pt();
    eleEta[nElectrons] = ele.eta();
    elePhi[nElectrons] = ele.phi();
    eleCharge[nElectrons] = ele.charge();
    eleE_SC[nElectrons] = 0;//ele.superCluster()->energy();
    eleEta_SC[nElectrons] = 0;//ele.superCluster()->eta();
    elePhi_SC[nElectrons] = 0;//ele.superCluster()->phi();
    eleSigmaIetaIeta[nElectrons] = ele.sigmaIetaIeta();
    eleFull5x5SigmaIetaIeta[nElectrons] = ele.full5x5_sigmaIetaIeta();
    eleR9[nElectrons] = ele.r9();
    ele_dEta[nElectrons] = ele.deltaEtaSuperClusterTrackAtVtx();
    ele_dPhi[nElectrons] = ele.deltaPhiSuperClusterTrackAtVtx();
    ele_HoverE[nElectrons] = ele.hcalOverEcal();
    ele_d0[nElectrons] = -ele.gsfTrack().get()->dxy(myPV->position());
    ele_dZ[nElectrons] = ele.gsfTrack().get()->dz(myPV->position());
    //ele_ip3d[nElectrons] = ele.dB(reco::Electron::PV3D);
    //ele_ip3dSignificance[nElectrons] = ele.dB(pat::Electron::PV3D)/ele.edB(pat::Electron::PV3D);
    ele_pileupIso[nElectrons] = ele.pfIsolationVariables().sumPUPt;
    ele_chargedIso[nElectrons] = ele.pfIsolationVariables().sumChargedHadronPt;
    ele_photonIso[nElectrons] = ele.pfIsolationVariables().sumPhotonEt;
    ele_neutralHadIso[nElectrons] = ele.pfIsolationVariables().sumNeutralHadronEt;
    //ele_MissHits[nElectrons] = ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

    //*************************************************                                                                                                 
    //Conversion Veto                                                                                                                                   
    //*************************************************                                                                                                 
    ele_PassConvVeto[nElectrons] = false;
    //if( beamSpot.isValid() && conversions.isValid() ) {
    //  ele_PassConvVeto[nElectrons] = !ConversionTools::hasMatchedConversion(ele,conversions,
    //                                                                        beamSpot->position());
    //} else {
    //  cout << "\n\nERROR!!! conversions not found!!!\n";
    //}

    // 1/E - 1/P                                                                                                                                        
    if( ele.ecalEnergy() == 0 ){
      ele_OneOverEminusOneOverP[nElectrons] = 1e30;
    } else if( !std::isfinite(ele.ecalEnergy())){
      ele_OneOverEminusOneOverP[nElectrons] = 1e30;
    } else {
      ele_OneOverEminusOneOverP[nElectrons] = 1./ele.ecalEnergy()  -  ele.eSuperClusterOverP()/ele.ecalEnergy();
    }

    //*************************************************                                                                                                 
    //ID MVA                                                                                                                                            
    ////*************************************************                                                                                                 
    //ele_IDMVATrig[nElectrons]    = myMVATrig->mvaValue(ele, false);
    //ele_IDMVANonTrig[nElectrons] = myMVANonTrig->mvaValue(ele,conversions, beamSpot->position(),false);
    //
    //ele_RegressionE[nElectrons] = ele.ecalRegressionEnergy();
    //ele_CombineP4[nElectrons]   = ele.ecalTrackRegressionEnergy();
    //
    //ele_ptrel[nElectrons]   = getLeptonPtRel( jets, &ele );
    //tuple<double,double,double> PFMiniIso = getPFMiniIsolation(packedPFCands, dynamic_cast<const reco::Candidate *>(&ele), 0.05, 0.2, 10., false, false);
    //ele_chargedMiniIso[nElectrons] = std::get<0>(PFMiniIso);
    //ele_photonAndNeutralHadronMiniIso[nElectrons] = std::get<1>(PFMiniIso);
    //ele_chargedPileupMiniIso[nElectrons] = std::get<2>(PFMiniIso);
    //ele_activityMiniIsoAnnulus[nElectrons] = ActivityPFMiniIsolationAnnulus( packedPFCands, dynamic_cast<const reco::Candidate *>(&ele), 0.4, 0.05, 0.2, 10.);
    //ele_passSingleEleTagFilter[nElectrons] = passSingleEleTagFilter;
    //ele_passTPOneTagFilter[nElectrons] = passTPOneTagFilter;
    //ele_passTPTwoTagFilter[nElectrons] = passTPTwoTagFilter;
    //ele_passTPOneProbeFilter[nElectrons] = passTPOneProbeFilter;
    //ele_passTPTwoProbeFilter[nElectrons] = passTPTwoProbeFilter;

    //missing some other stuff too
    
    nElectrons++;
  }

  return true;
}

bool RazorTuplizer::fillPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //noZS::EcalClusterLazyTools *lazyToolnoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, ebRecHitsToken_, eeRecHitsToken_);

  for (const reco::Photon &pho : *photons) {
    if (pho.pt() < 20) continue;

    //std::vector<float> vCov = lazyToolnoZS->localCovariances( *(pho.superCluster()->seed()) );

    //phoE[nPhotons] = pho.getCorrectedEnergy(reco::Photon::P4type::ecal_standard);
    phoPt[nPhotons] = pho.pt();
    //phoEta[nPhotons] = pho.eta(); //correct this for the vertex                                                                                       
    //phoPhi[nPhotons] = pho.phi(); //correct this for the vertex                                                                                       

    //phoSigmaIetaIeta[nPhotons] = pho.see();
    phoFull5x5SigmaIetaIeta[nPhotons] = pho.full5x5_sigmaIetaIeta();

    //phoR9[nPhotons] = pho.r9();                                                                                                                       
    //Use the noZS version of this according to Emanuele                                                                                                
    phoR9[nPhotons] = pho.full5x5_r9();

    pho_HoverE[nPhotons] = pho.hadTowOverEm();
    //pho_isConversion[nPhotons] = pho.hasConversionTracks();
    //pho_passEleVeto[nPhotons] = !hasMatchedPromptElectron(pho.superCluster(),electrons,conversions, beamSpot->position());

    //Don't use default miniAOD quantities for now                                                                                                      
    // pho_sumChargedHadronPt[nPhotons] = pho.chargedHadronIso();                                                                                       
    // pho_sumNeutralHadronEt[nPhotons] = pho.neutralHadronIso();                                                                                       
    // pho_sumPhotonEt[nPhotons] = pho.photonIso();                                                                                                     

    //**********************************************************                                                                                        
    //Compute PF isolation                                                                                                                              
    //absolute uncorrected isolations with footprint removal                                                                                            
    //**********************************************************                                                                                        
    //const float coneSizeDR = 0.3;
    //const float dxyMax = 0.1;
    //const float dzMax = 0.2;
    //float chargedIsoSum = 0;
    //float neutralHadronIsoSum = 0;
    //float photonIsoSum = 0;
    //
    //// First, find photon direction with respect to the good PV                                                                                         
    //math::XYZVector photon_directionWrtVtx(pho.superCluster()->x() - myPV->x(),
    //                                       pho.superCluster()->y() - myPV->y(),
    //                                       pho.superCluster()->z() - myPV->z());
    //// Loop over all PF candidates                                                                                                                      
    //for (const pat::PackedCandidate &candidate : *packedPFCands) {
    //
    //  // Check if this candidate is within the isolation cone                                                                                            
    //  float dR=deltaR(photon_directionWrtVtx.Eta(),photon_directionWrtVtx.Phi(),
    //                  candidate.eta(), candidate.phi());
    //  if( dR > coneSizeDR ) continue;
    //
    //  // Check if this candidate is not in the footprint                                                                                                
    //  bool inFootprint = false;
    //  for (auto itr : pho.associatedPackedPFCandidates()) {
    //    if ( &(*itr) == &candidate) {
    //      inFootprint = true;
    //    }
    //  }
    //  if( inFootprint ) continue;
    //
    //
    //  // Find candidate type                                                                                                                            
    //  reco::PFCandidate::ParticleType thisCandidateType = reco::PFCandidate::X;
    //
    //  // the neutral hadrons and charged hadrons can be of pdgId types                                                                                  
    //  // only 130 (K0L) and +-211 (pi+-) in packed candidates                                                                                           
    //  const int pdgId = candidate.pdgId();
    //  if( pdgId == 22 )
    //    thisCandidateType = reco::PFCandidate::gamma;
    //  else if( abs(pdgId) == 130) // PDG ID for K0L                                                                                                     
    //    thisCandidateType = reco::PFCandidate::h0;
    //  else if( abs(pdgId) == 211) // PDG ID for pi+-                                                                                                    
    //    thisCandidateType = reco::PFCandidate::h;
    //
    //
    //  // Increment the appropriate isolation sum                                                                                                        
    //  if( thisCandidateType == reco::PFCandidate::h ){
    //    // for charged hadrons, additionally check consistency                                                                                          
    //    // with the PV                                                                                                                                  
    //    float dxy = -999, dz = -999;
    //    dz = candidate.pseudoTrack().dz(myPV->position());
    //    dxy =candidate.pseudoTrack().dxy(myPV->position());
    //    if (fabs(dz) > dzMax) continue;
    //    if(fabs(dxy) > dxyMax) continue;
    //    // The candidate is eligible, increment the isolaiton                                                                                           
    //    chargedIsoSum += candidate.pt();
    //  }
    //  if( thisCandidateType == reco::PFCandidate::h0 )
    //    neutralHadronIsoSum += candidate.pt();
    //  if( thisCandidateType == reco::PFCandidate::gamma )
    //    photonIsoSum += candidate.pt();
    //}
    //pho_sumChargedHadronPt[nPhotons] = chargedIsoSum;
    //pho_sumNeutralHadronEt[nPhotons] = neutralHadronIsoSum;
    //pho_sumPhotonEt[nPhotons] = photonIsoSum;

    //*****************************************************************                                                                                 
    //Compute Worst Isolation Looping over all vertices                                                                                                 
    //*****************************************************************                                                                                 
    //const double ptMin = 0.0;
    //const float dRvetoBarrel = 0.0;
    //const float dRvetoEndcap = 0.0;
    //float dRveto = 0;
    //if (pho.isEB()) dRveto = dRvetoBarrel;
    //else dRveto = dRvetoEndcap;
    //
    //float worstIsolation = 999;
    //std::vector<float> allIsolations;
    //for(unsigned int ivtx=0; ivtx<vertices->size(); ++ivtx) {
    //
    //  // Shift the photon according to the vertex                                                                                                       
    //  reco::VertexRef vtx(vertices, ivtx);
    //  math::XYZVector photon_directionWrtVtx(pho.superCluster()->x() - vtx->x(),
    //                                         pho.superCluster()->y() - vtx->y(),
    //                                         pho.superCluster()->z() - vtx->z());
    //
    //  float sum = 0;
    //  // Loop over all PF candidates                                                                                                                    
    //  for (const pat::PackedCandidate &candidate : *packedPFCands) {
    //
    //    //require that PFCandidate is a charged hadron                                                                                                  
    //    const int pdgId = candidate.pdgId();
    //    if( abs(pdgId) != 211) continue;
    //
    //    if (candidate.pt() < ptMin)
    //      continue;
    //
    //    float dxy = -999, dz = -999;
    //    dz = candidate.pseudoTrack().dz(myPV->position());
    //    dxy =candidate.pseudoTrack().dxy(myPV->position());
    //    if( fabs(dxy) > dxyMax) continue;
    //    if ( fabs(dz) > dzMax) continue;
    //
    //    float dR = deltaR(photon_directionWrtVtx.Eta(), photon_directionWrtVtx.Phi(),
    //                      candidate.eta(),      candidate.phi());
    //    if(dR > coneSizeDR || dR < dRveto) continue;
    //
    //    sum += candidate.pt();
    //  }
    //
    //  allIsolations.push_back(sum);
    //}
    //
    //if( allIsolations.size()>0 )
    //  worstIsolation = * std::max_element( allIsolations.begin(), allIsolations.end() );
    //
    //pho_sumWorstVertexChargedHadronPt[nPhotons] = worstIsolation;

    //*****************************************************************                                                                                 
    //Photon ID MVA variable                                                                                                                            
    //*****************************************************************                                                                                 
    //pho_IDMVA[nPhotons] = myPhotonMVA->mvaValue( pho,  *rhoAll, photonIsoSum, chargedIsoSum, worstIsolation,
    //                                             lazyToolnoZS, false);
    //
    //pho_RegressionE[nPhotons] = pho.getCorrectedEnergy(reco::Photon::P4type::regression1);
    //pho_RegressionEUncertainty[nPhotons] = pho.getCorrectedEnergyError(reco::Photon::P4type::regression1);
    //
    ////compute photon corrected 4-mometum                                                                                                                
    //TVector3 phoPos( pho.superCluster()->x(), pho.superCluster()->y(), pho.superCluster()->z() );
    //TVector3 vtxPos( pvX, pvY, pvZ );
    //TLorentzVector phoP4 = photonP4FromVtx( vtxPos, phoPos, pho_RegressionE[nPhotons] );
    //phoEta[nPhotons] = phoP4.Eta();
    //phoPhi[nPhotons] = phoP4.Phi();
    //
    //pho_superClusterEta[nPhotons] = pho.superCluster()->eta();
    //pho_superClusterPhi[nPhotons] = pho.superCluster()->phi();
    //pho_hasPixelSeed[nPhotons] = pho.hasPixelSeed();

    nPhotons++;
  }

  //delete lazyToolnoZS;
  return true;
}

bool RazorTuplizer::fillJets(){
  for (const reco::PFJet &j : *jets) {
    // UNCORRECTED JETS!
    if (j.pt() < 20) continue;
    jetE[nJets] = j.energy();
    jetPt[nJets] = j.pt();
    jetEta[nJets] = j.eta();
    jetPhi[nJets] = j.phi();
    //jetCSV[nJets] = j.bDiscriminator("pfCombinedSecondaryVertexBJetTags");
    //jetCISV[nJets] = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    jetMass[nJets] = j.mass();
    //jetJetArea[nJets] = j.jetArea();
    //jetPileupE[nJets] = j.pileup();
    //jetPileupId[nJets] = j.userFloat("pileupJetId:fullDiscriminant");
    jetPileupIdFlag[nJets] = 0;
    jetPassIDLoose[nJets] = true;
    jetPassIDTight[nJets] = true;
    jetPassMuFrac[nJets]  = ( j.muonEnergyFraction() < 0.80 );
    jetPassEleFrac[nJets]  = ( j.electronEnergyFraction() < 0.90 );
    //if (useGen_) {
    //  jetPartonFlavor[nJets] = j.partonFlavour();
    //  jetHadronFlavor[nJets] = j.hadronFlavour();
    //} else {
    jetPartonFlavor[nJets] = -999;
    jetHadronFlavor[nJets] = -999;
    //}

    //extra jet information (may be only needed for debugging)
    //jetChargedEMEnergyFraction[nJets] = j.chargedEmEnergyFraction();
    //jetNeutralEMEnergyFraction[nJets] = j.neutralEmEnergyFraction();
    //jetChargedHadronEnergyFraction[nJets] = j.chargedHadronEnergyFraction();
    //jetNeutralHadronEnergyFraction[nJets] = j.neutralHadronEnergyFraction();
    //jetMuonEnergyFraction[nJets] =  j.muonEnergyFraction();
    //jetHOEnergyFraction[nJets] =  j.hoEnergyFraction();
    //jetHFHadronEnergyFraction[nJets] =  j.HFHadronEnergyFraction();
    //jetHFEMEnergyFraction[nJets] =  j.HFEMEnergyFraction();

    nJets++;
  }

  return true;
}

bool RazorTuplizer::fillJetsAK8(){
  // UNCORRECTED
  for (const reco::PFJet &j : *jetsAK8) {
    fatJetE[nFatJets] = j.energy();
    fatJetPt[nFatJets] = j.pt();
    fatJetEta[nFatJets] = j.eta();
    fatJetPhi[nFatJets] = j.phi();
    //fatJetPrunedM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSPrunedLinks");
    //fatJetTrimmedM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSTrimmedLinks");
    //fatJetFilteredM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSFilteredLinks");
    //fatJetTau1[nFatJets] =  (float) j.userFloat("NjettinessAK8:tau1");
    //fatJetTau2[nFatJets] =  (float) j.userFloat("NjettinessAK8:tau2");
    //fatJetTau3[nFatJets] =  (float) j.userFloat("NjettinessAK8:tau3");
    nFatJets++;
  }

  return true;
}

bool RazorTuplizer::fillMC(){

  for(const reco::GenJet &j : *genJets){
    genJetE[nGenJets] = j.energy();
    genJetPt[nGenJets] = j.pt();
    genJetEta[nGenJets] = j.eta();
    genJetPhi[nGenJets] = j.phi();
    nGenJets++;
  }
  
  //const reco::MET &Met = mets->front();
  //genMetPt = Met.genMET()->pt();
  //genMetPhi = Met.genMET()->phi();
  
  bool foundGenVertex = false;
  for(size_t i=0; i<genParticles->size();i++){
    if (!foundGenVertex) {
      for (unsigned int j=0; j<(*genParticles)[i].numberOfDaughters(); ++j) {
	const reco::Candidate *dau = (*genParticles)[i].daughter(j);
	if (dau) {
	  genVertexX = dau->vx();
	  genVertexY = dau->vy();
	  genVertexZ = dau->vz();
	  foundGenVertex = true;
	  break;
	}
      }
    }
  }
  
  genWeight = genInfo->weight();
  genSignalProcessID = genInfo->signalProcessID();
  genQScale = genInfo->qScale();
  genAlphaQCD = genInfo->alphaQCD();
  genAlphaQED = genInfo->alphaQED();

  return true;
}

bool RazorTuplizer::fillGenParticles(){
  std::vector<const reco::Candidate*> prunedV;//Allows easier comparison for mother finding                                                             
  //Fills selected gen particles                                                                                                                        
  for(size_t i=0; i<genParticles->size();i++){
    if(
       (abs((*genParticles)[i].pdgId()) >= 1 && abs((*genParticles)[i].pdgId()) <= 6
        && ( (*genParticles)[i].status() < 30
             )
        )
       || (abs((*genParticles)[i].pdgId()) >= 11 && abs((*genParticles)[i].pdgId()) <= 16)
       || (abs((*genParticles)[i].pdgId()) == 21
           && (*genParticles)[i].status() < 30
           )
       || (abs((*genParticles)[i].pdgId()) >= 22 && abs((*genParticles)[i].pdgId()) <= 25
           && ( (*genParticles)[i].status() < 30
                )
           )
       || (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)
       || (abs((*genParticles)[i].pdgId()) >= 1000001 && abs((*genParticles)[i].pdgId()) <= 1000039)
       ){
      prunedV.push_back(&(*genParticles)[i]);
    }
  }

  //Total number of gen particles                                                                                                                       
  nGenParticle = prunedV.size();
  //Look for mother particle and Fill gen variables                                                                                                     
  for(unsigned int i = 0; i < prunedV.size(); i++){
    gParticleId[i] = prunedV[i]->pdgId();
    gParticleStatus[i] = prunedV[i]->status();
    gParticleE[i] = prunedV[i]->energy();
    gParticlePt[i] = prunedV[i]->pt();
    gParticleEta[i] = prunedV[i]->eta();
    gParticlePhi[i] = prunedV[i]->phi();
    gParticleMotherId[i] = 0;
    gParticleMotherIndex[i] = -1;
    //if(prunedV[i]->numberOfMothers() > 0){

      //find the ID of the first mother that has a different ID than the particle itself                                                                
      //const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(prunedV[i]);
      //if (firstMotherWithDifferentID) {
      //gParticleMotherId[i] = firstMotherWithDifferentID->pdgId();
      //}

      //find the mother and keep going up the mother chain if the ID's are the same                                                                     
      //const reco::Candidate* originalMotherWithSameID = findOriginalMotherWithSameID(prunedV[i]);
      //for(unsigned int j = 0; j < prunedV.size(); j++){
      //if(prunedV[j] == originalMotherWithSameID){
      //gParticleMotherIndex[i] = j;
      //break;
      //}
      //}
    //} else {
    gParticleMotherIndex[i] = -1;
    //}
  }
  return true;
}

void RazorTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  //initialize                                                                                                                                      
  resetBranches();
  loadEvent(iEvent); //loads objects and resets tree branches

  NEvents->Fill(0);

  bool isGoodEvent = 
    fillEventInfo(iEvent) 
    && fillMuons()
    && fillElectrons()
    && fillPhotons(iEvent,iSetup)
    && fillJets()
    && fillJetsAK8();

  bool isGoodMCEvent = true;
  if (useGen_) {
    isGoodMCEvent = fillMC() 
      && fillPileUp() 
      && fillGenParticles();
  }
  isGoodEvent = isGoodEvent&&isGoodMCEvent;

  //fill the tree if the event wasn't rejected                                                                                                      
  if(isGoodEvent) RazorEvents->Fill();
}

//------ Method called once each job just before starting event loop ------//                                                                       
void RazorTuplizer::beginJob(){
  setBranches();
}

//------ Method called once each job just after ending the event loop ------//                                                                      
void RazorTuplizer::endJob(){
}

//define this as a plug-in                                                                                                                          
DEFINE_FWK_MODULE(RazorTuplizer);
