// babymaker: Makes flat ntuples with HT, MET, MC, jets, and leptons

#include <string>
#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "hltstudy/babymaker/interface/babymaker.h"

#include "TMath.h"

namespace{
  template<typename T>
    std::auto_ptr<T> make_auto(T * const p=nullptr){
    return std::auto_ptr<T>{p};
  }
}

babymaker::babymaker(const edm::ParameterSet& iConfig) {

  produces<float> ("wl1ht200").setBranchAlias("wl1ht200");
  produces<float> ("pfht").setBranchAlias("pf_ht");
  produces<float> ("caloht").setBranchAlias("calo_ht");
  produces<float> ("genht").setBranchAlias("gen_ht");

  produces<float> ("genmet").setBranchAlias("gen_met");
  produces<float> ("genmetcalo").setBranchAlias("gen_metcalo");
  produces<float> ("genmetcalononprompt").setBranchAlias("gen_metcalononprompt");

  produces<std::vector<float> > ("elept").setBranchAlias("ele_pt");
  produces<std::vector<float> > ("eleeta").setBranchAlias("ele_eta");
  produces<std::vector<float> > ("elephi").setBranchAlias("ele_phi");

  produces<std::vector<float> > ("elspt").setBranchAlias("els_pt");
  produces<std::vector<float> > ("elseta").setBranchAlias("els_eta");
  produces<std::vector<float> > ("elsphi").setBranchAlias("els_phi");
  produces<std::vector<float> > ("elsecaliso").setBranchAlias("els_ecal_iso");
  produces<std::vector<float> > ("elshcaliso").setBranchAlias("els_hcal_iso");
  produces<std::vector<float> > ("elstrackiso").setBranchAlias("els_track_iso");
    
  produces<std::vector<float> > ("muspt").setBranchAlias("mus_pt");
  produces<std::vector<float> > ("museta").setBranchAlias("mus_eta");
  produces<std::vector<float> > ("musphi").setBranchAlias("mus_phi");
  produces<std::vector<float> > ("musiso").setBranchAlias("mus_iso");
    
  produces<std::vector<float> > ("pfjetspt").setBranchAlias("pfjets_pt");
  produces<std::vector<float> > ("pfjetseta").setBranchAlias("pfjets_eta");
  produces<std::vector<float> > ("pfjetsphi").setBranchAlias("pfjets_phi");

  produces<std::vector<float> > ("genjetspt").setBranchAlias("genjets_pt");
  produces<std::vector<float> > ("genjetseta").setBranchAlias("genjets_eta");
  produces<std::vector<float> > ("genjetsphi").setBranchAlias("genjets_phi");
    
  produces<std::vector<float> > ("genelspt").setBranchAlias("genels_pt");
  produces<std::vector<float> > ("genelseta").setBranchAlias("genels_eta");
  produces<std::vector<float> > ("genelsphi").setBranchAlias("genels_phi");
  produces<std::vector<int> >   ("genelsmomid").setBranchAlias("genels_mom_id");
  produces<std::vector<int> >   ("genelsgmomid").setBranchAlias("genels_gmom_id");
  produces<std::vector<int> >   ("genelsggmomid").setBranchAlias("genels_ggmom_id");
    
  produces<std::vector<float> > ("genmuspt").setBranchAlias("genmus_pt");
  produces<std::vector<float> > ("genmuseta").setBranchAlias("genmus_eta");
  produces<std::vector<float> > ("genmusphi").setBranchAlias("genmus_phi");
  produces<std::vector<int> >   ("genmusmomid").setBranchAlias("genmus_mom_id");
  produces<std::vector<int> >   ("genmusgmomid").setBranchAlias("genmus_gmom_id");
  produces<std::vector<int> >   ("genmusggmomid").setBranchAlias("genmus_ggmom_id");

  produces<std::vector<float> > ("recomuspt") .setBranchAlias("reco_mus_pt");
  produces<std::vector<float> > ("recomuseta").setBranchAlias("reco_mus_eta");
  produces<std::vector<float> > ("recomusphi").setBranchAlias("reco_mus_phi");

  produces<std::vector<float> > ("recoelspt") .setBranchAlias("reco_els_pt");
  produces<std::vector<float> > ("recoelseta").setBranchAlias("reco_els_eta");
  produces<std::vector<float> > ("recoelsphi").setBranchAlias("reco_els_phi");

  produces<std::vector<float> > ("recojetpt") .setBranchAlias("reco_jet_pt");
  produces<std::vector<float> > ("recojeteta").setBranchAlias("reco_jet_eta");
  produces<std::vector<float> > ("recojetphi").setBranchAlias("reco_jet_phi");

  produces<std::vector<float> > ("recogenjetpt") .setBranchAlias("reco_genjet_pt");
  produces<std::vector<float> > ("recogenjeteta").setBranchAlias("reco_genjet_eta");
  produces<std::vector<float> > ("recogenjetphi").setBranchAlias("reco_genjet_phi");
  
  produces<std::vector<float> > ("bjetspt").setBranchAlias("bjets_pt");
  produces<std::vector<float> > ("bjetseta").setBranchAlias("bjets_eta");
  produces<std::vector<float> > ("bjetsphi").setBranchAlias("bjets_phi");
  produces<std::vector<float> > ("bjetscsv").setBranchAlias("bjets_csv");


  produces<float> ("metpt").setBranchAlias("met_pt");
  produces<float> ("meteta").setBranchAlias("met_eta");
  produces<float> ("metphi").setBranchAlias("met_phi");

  produces<float> ("pfmhtpt").setBranchAlias("pf_mht_pt");
  produces<float> ("pfmhteta").setBranchAlias("pf_mht_eta");
  produces<float> ("pfmhtphi").setBranchAlias("pf_mht_phi");

  produces<float> ("calomhtpt").setBranchAlias("calo_mht_pt");
  produces<float> ("calomhteta").setBranchAlias("calo_mht_eta");
  produces<float> ("calomhtphi").setBranchAlias("calo_mht_phi");

  produces<unsigned int> ("ngenlep").setBranchAlias("ngenlep");

  hltPfJetsInputTag = iConfig.getParameter<edm::InputTag>("hltPfJetsInputTag_");
  hltElectronInputString = iConfig.getParameter<std::string>("hltElectronInputString_");
  hltElectronPtInputTag = edm::InputTag(hltElectronInputString, "elept");
  hltElectronPhiInputTag = edm::InputTag(hltElectronInputString, "elephi");
  hltElectronEtaInputTag = edm::InputTag(hltElectronInputString, "eleeta");
  hltElectronEcalIsoInputTag = edm::InputTag(hltElectronInputString, "eleecaliso");
  hltElectronHcalIsoInputTag = edm::InputTag(hltElectronInputString, "elehcaliso");
  hltElectronTrackIsoInputTag = edm::InputTag(hltElectronInputString, "eletrackiso");
  hltMuonInputString = iConfig.getParameter<std::string>("hltMuonInputString_");
  hltMuonPtInputTag = edm::InputTag(hltMuonInputString, "mupt");
  hltMuonPhiInputTag = edm::InputTag(hltMuonInputString, "muphi");
  hltMuonEtaInputTag = edm::InputTag(hltMuonInputString, "mueta");
  hltMuonIsoInputTag = edm::InputTag(hltMuonInputString, "muiso");
  hltPfMetInputTag = iConfig.getParameter<edm::InputTag>("hltPfMetInputTag_");
  hltPfHTInputTag = iConfig.getParameter<edm::InputTag>("hltPfHTInputTag_");
  hltCaloHTInputTag = iConfig.getParameter<edm::InputTag>("hltCaloHTInputTag_");
  hltGenJetsInputTag = iConfig.getParameter<edm::InputTag>("hltGenJetsInputTag_");

  recoPfJetsInputTag = iConfig.getParameter<edm::InputTag>("recoPfJetsInputTag_");
  recoElectronInputTag = iConfig.getParameter<edm::InputTag>("recoElectronInputTag_");
  recoMuonInputTag = iConfig.getParameter<edm::InputTag>("recoMuonInputTag_");
  recoPfMetInputTag = iConfig.getParameter<edm::InputTag>("recoPfMetInputTag_");
  recoGenJetsInputTag = iConfig.getParameter<edm::InputTag>("recoGenJetsInputTag_");
  
  // btagging
  m_Jets   =  iConfig.getParameter<edm::InputTag>("m_Jets_");
  m_JetTags = iConfig.getParameter<edm::InputTag>("m_JetTags_"); 
  m_JetsToken = consumes<std::vector<reco::CaloJet> >(m_Jets);
  m_JetTagsToken = consumes<reco::JetTagCollection>(m_JetTags);
}

babymaker::~babymaker() {}

void  babymaker::beginJob() {
}

void babymaker::endJob() {
}

// ------------ method called to produce the data  ------------
void babymaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto els_pt  = make_auto(new std::vector<float>);
  auto els_eta = make_auto(new std::vector<float>);
  auto els_phi = make_auto(new std::vector<float>);
  auto els_ecal_iso = make_auto(new std::vector<float>);
  auto els_hcal_iso = make_auto(new std::vector<float>);
  auto els_track_iso = make_auto(new std::vector<float>);

  auto ele_pt = make_auto(new std::vector<float>);
  auto ele_phi = make_auto(new std::vector<float>);
  auto ele_eta = make_auto(new std::vector<float>);
  
  auto mus_pt  = make_auto(new std::vector<float>);
  auto mus_eta = make_auto(new std::vector<float>);
  auto mus_phi = make_auto(new std::vector<float>);
  auto mus_iso = make_auto(new std::vector<float>);
    
  auto pfjets_pt  = make_auto(new std::vector<float>);
  auto pfjets_eta = make_auto(new std::vector<float>);
  auto pfjets_phi = make_auto(new std::vector<float>);

  auto genjets_pt  = make_auto(new std::vector<float>);
  auto genjets_eta = make_auto(new std::vector<float>);
  auto genjets_phi = make_auto(new std::vector<float>);

  auto genels_pt       = make_auto(new std::vector<float>);
  auto genels_eta      = make_auto(new std::vector<float>);
  auto genels_phi      = make_auto(new std::vector<float>);
  auto genels_mom_id   = make_auto(new std::vector<int>);
  auto genels_gmom_id  = make_auto(new std::vector<int>);
  auto genels_ggmom_id = make_auto(new std::vector<int>);

  auto genmus_pt       = make_auto(new std::vector<float>);
  auto genmus_eta      = make_auto(new std::vector<float>);
  auto genmus_phi      = make_auto(new std::vector<float>);
  auto genmus_mom_id   = make_auto(new std::vector<int>);
  auto genmus_gmom_id  = make_auto(new std::vector<int>);
  auto genmus_ggmom_id = make_auto(new std::vector<int>);

  auto reco_mus_pt  = make_auto(new std::vector<float>);
  auto reco_mus_eta = make_auto(new std::vector<float>);
  auto reco_mus_phi = make_auto(new std::vector<float>);

  auto reco_els_pt  = make_auto(new std::vector<float>);
  auto reco_els_eta = make_auto(new std::vector<float>);
  auto reco_els_phi = make_auto(new std::vector<float>);

  auto reco_jet_pt  = make_auto(new std::vector<float>);
  auto reco_jet_eta = make_auto(new std::vector<float>);
  auto reco_jet_phi = make_auto(new std::vector<float>);

  auto reco_genjet_pt  = make_auto(new std::vector<float>);
  auto reco_genjet_eta = make_auto(new std::vector<float>);
  auto reco_genjet_phi = make_auto(new std::vector<float>);
  
  auto bjets_pt     = make_auto(new std::vector<float>);
  auto bjets_eta    = make_auto(new std::vector<float>);
  auto bjets_phi    = make_auto(new std::vector<float>);
  auto bjets_csv    = make_auto(new std::vector<float>);

  auto met_pt  = make_auto(new float);
  auto met_eta = make_auto(new float);
  auto met_phi = make_auto(new float);

  auto pf_ht    = make_auto(new float);
  auto gen_ht   = make_auto(new float);
  auto wl1ht200 = make_auto(new float);

  auto gen_met              = make_auto(new float);
  auto gen_metcalo          = make_auto(new float);
  auto gen_metcalononprompt = make_auto(new float);

  auto pf_mht_pt  = make_auto(new float);
  auto pf_mht_eta = make_auto(new float);
  auto pf_mht_phi = make_auto(new float);

  auto calo_ht = make_auto(new float);
  auto calo_mht_pt = make_auto(new float);
  auto calo_mht_eta = make_auto(new float);
  auto calo_mht_phi = make_auto(new float);

  auto pseudoMT2_hlt = make_auto(new float);
  auto pseudoMT2_snt = make_auto(new float);
  auto mt2_hlt       = make_auto(new float);
  auto mt2_snt       = make_auto(new float);

  auto ngenlep = make_auto(new unsigned int);

  edm::Handle<edm::View<reco::PFJet> > jet_h;
  if(hltPfJetsInputTag.label() != "unused"){
    iEvent.getByLabel(hltPfJetsInputTag, jet_h);
    for(auto jet_it = jet_h->begin(); jet_it != jet_h->end(); ++jet_it){
      if(jet_it->pt() < 35.0) continue;
      pfjets_pt  ->push_back(jet_it->pt());
      pfjets_eta ->push_back(jet_it->eta());
      pfjets_phi ->push_back(jet_it->phi());
    } 
  }

  edm::Handle<edm::View<reco::MET> > pfht_h;
  if(hltPfHTInputTag.label() != "unused"){
    iEvent.getByLabel(hltPfHTInputTag, pfht_h);
    *pf_ht    = (pfht_h->front()).sumEt();
    *pf_mht_pt    = (pfht_h->front()).pt();
    *pf_mht_eta   = (pfht_h->front()).eta();
    *pf_mht_phi   = (pfht_h->front()).phi();
  }

  edm::Handle<edm::View<reco::MET> > met_h;
  if(hltPfMetInputTag.label() != "unused"){
    iEvent.getByLabel(hltPfMetInputTag, met_h);
    *met_pt   = (met_h->front()).pt();
    *met_eta  = (met_h->front()).eta();
    *met_phi  = (met_h->front()).phi();
  }

  if(hltCaloHTInputTag.label() != "unused"){
    edm::Handle<edm::View<reco::MET> > calo_ht_h;
    iEvent.getByLabel(hltCaloHTInputTag, calo_ht_h);
    *calo_ht    = (calo_ht_h->front()).sumEt();
    *calo_mht_pt    = (calo_ht_h->front()).pt();
    *calo_mht_eta   = (calo_ht_h->front()).eta();
    *calo_mht_phi   = (calo_ht_h->front()).phi();
  }

  edm::Handle<edm::View<reco::GenJet> > genjet_h;
  iEvent.getByLabel(hltGenJetsInputTag, genjet_h);

  edm::Handle<reco::GenParticleCollection> genparts;
  iEvent.getByLabel("genParticles", genparts);


  edm::Handle<edm::View<reco::GenMET> >genmet_h;
  iEvent.getByLabel("genMetTrue", genmet_h);
  
  edm::Handle<edm::View<reco::GenMET> >genmetcalo_h;
  iEvent.getByLabel("genMetCalo", genmetcalo_h);
  
  edm::Handle<edm::View<reco::GenMET> >genmetcalononprompt_h;
  iEvent.getByLabel("genMetCaloAndNonPrompt", genmetcalononprompt_h);

  //edm::Handle<reco::JetTagCollection> btag_h ;
  //iEvent.getByLabel("hltL3CombinedSecondaryVertexBJetTags", btag_h) ;
  
  edm::Handle<reco::MuonCollection> reco_muon_h;
  if(recoMuonInputTag.label() != "unused"){
    iEvent.getByLabel(recoMuonInputTag, reco_muon_h);
    for(auto reco_muon_it = reco_muon_h->begin(); reco_muon_it!=reco_muon_h->end(); ++reco_muon_it){
      reco_mus_pt ->push_back(reco_muon_it->pt());
      reco_mus_eta->push_back(reco_muon_it->eta());
      reco_mus_phi->push_back(reco_muon_it->phi());
    }
  }

  edm::Handle<reco::GsfElectronCollection> reco_electron_h;
  if(recoElectronInputTag.label() != "unused"){
    iEvent.getByLabel(recoElectronInputTag, reco_electron_h);
    for(auto reco_electron_it = reco_electron_h->begin(); reco_electron_it!=reco_electron_h->end(); ++reco_electron_it){
      reco_els_pt ->push_back(reco_electron_it->pt());
      reco_els_eta->push_back(reco_electron_it->eta());
      reco_els_phi->push_back(reco_electron_it->phi());
    }
  }

  edm::Handle<reco::PFJetCollection> reco_jet_h;
  if(recoPfJetsInputTag.label() != "unused"){
    iEvent.getByLabel(recoPfJetsInputTag, reco_jet_h);
    for(auto reco_jet_it = reco_jet_h->begin(); reco_jet_it!=reco_jet_h->end(); ++reco_jet_it){
      if(reco_jet_it->pt()<35.0) continue;
      reco_jet_pt->push_back(reco_jet_it->pt());
      reco_jet_eta->push_back(reco_jet_it->eta());
      reco_jet_phi->push_back(reco_jet_it->phi());
    }
  }

  edm::Handle<reco::GenJetCollection> reco_genjet_h;
  if(recoGenJetsInputTag.label() != "unused"){
    iEvent.getByLabel(recoGenJetsInputTag, reco_genjet_h);
    for(auto reco_genjet_it = reco_genjet_h->begin(); reco_genjet_it!=reco_genjet_h->end(); ++reco_genjet_it){
      if(reco_genjet_it->pt()<35.0) continue;
      reco_genjet_pt->push_back(reco_genjet_it->pt());
      reco_genjet_eta->push_back(reco_genjet_it->eta());
      reco_genjet_phi->push_back(reco_genjet_it->phi());
    }
  }

  *ngenlep = 0;
  for (auto genp = genparts->cbegin(); genp != genparts->cend(); ++genp ) {
    if ((genp->status() == 3 || genp->status() == 23) && (abs(genp->pdgId()) == 11 || abs(genp->pdgId()) == 13 || abs(genp->pdgId()) == 15)) (*ngenlep)++;
        
    //if (genp->status() != 3 && genp->status() != 23) continue;  

    int id = genp->pdgId();
    const reco::Candidate* mom = genp->mother();
    int mom_id = mom ? mom->pdgId() : 0;
    const reco::Candidate* gmom = mom ? mom->mother() : nullptr;
    int gmom_id = gmom ? gmom->pdgId() : 0;
    const reco::Candidate* ggmom = gmom ? gmom->mother() : nullptr;
    int ggmom_id = ggmom ? ggmom->pdgId() : 0;

        
    if( TMath::Abs(id) != 11 && TMath::Abs(id) != 13 ) continue; 
    if( TMath::Abs(mom_id) != 24 && TMath::Abs(gmom_id) != 24 && TMath::Abs(ggmom_id) != 24 
	&& genp->status() != 3 && genp->status() != 23 ) continue;
        
    if( TMath::Abs(id) == 11 ) {
      genels_pt  ->push_back(genp->p4().pt());
      genels_eta ->push_back(genp->p4().eta());
      genels_phi ->push_back(genp->p4().phi());
      genels_mom_id->push_back(mom_id);
      genels_gmom_id->push_back(gmom_id);
      genels_ggmom_id->push_back(ggmom_id);
    }
    if( TMath::Abs(id) == 13 ) {
      genmus_pt  ->push_back(genp->p4().pt());
      genmus_eta ->push_back(genp->p4().eta());
      genmus_phi ->push_back(genp->p4().phi());
      genmus_mom_id->push_back(mom_id);
      genmus_gmom_id->push_back(gmom_id);
      genmus_ggmom_id->push_back(ggmom_id);
    }
  }

  if(hltElectronInputString != "unused"){
    edm::Handle<std::vector<float> > els_pt_h;
    edm::Handle<std::vector<float> > els_phi_h;
    edm::Handle<std::vector<float> > els_eta_h;
    edm::Handle<std::vector<float> > els_ecal_iso_h;
    edm::Handle<std::vector<float> > els_hcal_iso_h;
    edm::Handle<std::vector<float> > els_track_iso_h;
    iEvent.getByLabel(hltElectronPtInputTag, els_pt_h);
    iEvent.getByLabel(hltElectronPhiInputTag, els_phi_h);
    iEvent.getByLabel(hltElectronEtaInputTag, els_eta_h);
    iEvent.getByLabel(hltElectronEcalIsoInputTag, els_ecal_iso_h);
    iEvent.getByLabel(hltElectronHcalIsoInputTag, els_hcal_iso_h);
    iEvent.getByLabel(hltElectronTrackIsoInputTag, els_track_iso_h);
    for(const auto& pt: *els_pt_h) els_pt->push_back(pt);
    for(const auto& phi: *els_phi_h) els_phi->push_back(phi);
    for(const auto& eta: *els_eta_h) els_eta->push_back(eta);
    for(const auto& iso: *els_ecal_iso_h) els_ecal_iso->push_back(iso);
    for(const auto& iso: *els_hcal_iso_h) els_hcal_iso->push_back(iso);
    for(const auto& iso: *els_track_iso_h) els_track_iso->push_back(iso);
  }

  if(hltMuonInputString !="unused") {
    edm::Handle<std::vector<float> > mus_pt_h;
    edm::Handle<std::vector<float> > mus_phi_h;
    edm::Handle<std::vector<float> > mus_eta_h;
    edm::Handle<std::vector<float> > mus_iso_h;
    iEvent.getByLabel(hltMuonPtInputTag, mus_pt_h);
    iEvent.getByLabel(hltMuonPhiInputTag, mus_phi_h);
    iEvent.getByLabel(hltMuonEtaInputTag, mus_eta_h);
    iEvent.getByLabel(hltMuonIsoInputTag, mus_iso_h);
    for(const auto& pt: *mus_pt_h) mus_pt->push_back(pt);
    for(const auto& phi: *mus_phi_h) mus_phi->push_back(phi);
    for(const auto& eta: *mus_eta_h) mus_eta->push_back(eta);
    for(const auto& iso: *mus_iso_h) mus_iso->push_back(iso);
  }

  *gen_ht = 0;
  for(auto genjet_it = genjet_h->begin(); genjet_it != genjet_h->end(); ++genjet_it){

    if(genjet_it->pt() < 35.0) continue;

    genjets_pt  ->push_back(genjet_it->pt());
    genjets_eta ->push_back(genjet_it->eta());
    genjets_phi ->push_back(genjet_it->phi());

    if(genjet_it->pt()>40 && genjet_it->eta() < 3) *gen_ht += genjet_it->pt();

  }

 
  reco::GenMET genmet_obj(genmet_h->at(0));
  *gen_met = genmet_obj.pt();
  
  reco::GenMET genmetcalo_obj(genmetcalo_h->at(0));
  *gen_metcalo = genmetcalo_obj.pt();
  
  reco::GenMET genmetcalononprompt_obj(genmetcalononprompt_h->at(0));
  *gen_metcalononprompt = genmetcalononprompt_obj.pt();

  // Btagging --------------------------------------- begin
  if(m_Jets.label() != "unused" && m_JetTags.label() != "unused"){
    typedef std::vector<reco::CaloJet> TCollection;
    typedef edm::Ref<TCollection> TRef;
    
    edm::Handle<TCollection> h_Jets;
    iEvent.getByToken(m_JetsToken, h_Jets);
    
    edm::Handle<reco::JetTagCollection> h_JetTags;
    iEvent.getByToken(m_JetTagsToken, h_JetTags);
    
    TRef jetRef;
    
    // Look at all jets in decreasing order of Et.
    for (reco::JetTagCollection::const_iterator jet = h_JetTags->begin(); jet != h_JetTags->end(); ++jet) {
      
      jetRef = TRef(h_Jets,jet->first.key());
      if(jet->first->pt() < 35) continue;
      
      if(0) 
	{
	  std::cout << "Jet pT = " << jet->first->pt()
		    << ", CSV = " << jet->second << std::endl;
	}
      
      bjets_pt ->push_back( jet->first->pt() );
      bjets_eta ->push_back( jet->first->eta() );
      bjets_phi ->push_back( jet->first->phi() );
      bjets_csv ->push_back( jet->second );
      
    }
  }else{
    bjets_pt->clear();
    bjets_eta->clear();
    bjets_phi->clear();
    bjets_csv->clear();
  }
 
  // Btagging --------------------------------------- end 

    



  // Turn on curve for L1 HTT200 (calculated by Dominick)
  *wl1ht200 = (0.5*TMath::Erf((1.35121e-02)*(*gen_ht-(3.02695e+02)))+0.5);



  iEvent.put(els_pt,   "elspt" );
  iEvent.put(els_eta,  "elseta" );
  iEvent.put(els_phi,  "elsphi" );
  iEvent.put(els_ecal_iso,  "elsecaliso" );
  iEvent.put(els_hcal_iso,  "elshcaliso" );
  iEvent.put(els_track_iso,  "elstrackiso" );

  iEvent.put(ele_pt,   "elept" );
  iEvent.put(ele_eta,  "eleeta" );
  iEvent.put(ele_phi,  "elephi" );
    
  iEvent.put(mus_pt,   "muspt" );
  iEvent.put(mus_eta,  "museta" );
  iEvent.put(mus_phi,  "musphi" );
  iEvent.put(mus_iso,  "musiso" );

  iEvent.put(pfjets_pt,   "pfjetspt" );
  iEvent.put(pfjets_eta,  "pfjetseta" );
  iEvent.put(pfjets_phi,  "pfjetsphi" );

  iEvent.put(genjets_pt,   "genjetspt" );
  iEvent.put(genjets_eta,  "genjetseta" );
  iEvent.put(genjets_phi,  "genjetsphi" );

  iEvent.put(genels_pt,       "genelspt" );
  iEvent.put(genels_eta,      "genelseta" );
  iEvent.put(genels_phi,      "genelsphi" );
  iEvent.put(genels_mom_id,   "genelsmomid");
  iEvent.put(genels_gmom_id,  "genelsgmomid");
  iEvent.put(genels_ggmom_id, "genelsggmomid");

  iEvent.put(genmus_pt,       "genmuspt" );
  iEvent.put(genmus_eta,      "genmuseta" );
  iEvent.put(genmus_phi,      "genmusphi" );
  iEvent.put(genmus_mom_id,   "genmusmomid");
  iEvent.put(genmus_gmom_id,  "genmusgmomid");
  iEvent.put(genmus_ggmom_id, "genmusggmomid");

  iEvent.put(reco_mus_pt,  "recomuspt");
  iEvent.put(reco_mus_eta, "recomuseta");
  iEvent.put(reco_mus_phi, "recomusphi");

  iEvent.put(reco_els_pt,  "recoelspt");
  iEvent.put(reco_els_eta, "recoelseta");
  iEvent.put(reco_els_phi, "recoelsphi");

  iEvent.put(reco_jet_pt,  "recojetpt");
  iEvent.put(reco_jet_eta, "recojeteta");
  iEvent.put(reco_jet_phi, "recojetphi");

  iEvent.put(reco_genjet_pt,  "recogenjetpt");
  iEvent.put(reco_genjet_eta, "recogenjeteta");
  iEvent.put(reco_genjet_phi, "recogenjetphi");
  
  iEvent.put(bjets_pt,      "bjetspt" ); 
  iEvent.put(bjets_eta,     "bjetseta" ); 
  iEvent.put(bjets_phi,     "bjetsphi" ); 
  iEvent.put(bjets_csv,     "bjetscsv" ); 

  iEvent.put(met_pt,   "metpt" );
  iEvent.put(met_eta,  "meteta" );
  iEvent.put(met_phi,  "metphi" );

  iEvent.put(pf_ht,   "pfht" );
  iEvent.put(gen_ht,   "genht" );
  iEvent.put(wl1ht200,   "wl1ht200" );

  iEvent.put(gen_met, "genmet");
  iEvent.put(gen_metcalo, "genmetcalo");
  iEvent.put(gen_metcalononprompt, "genmetcalononprompt");

  iEvent.put(pf_mht_pt,   "pfmhtpt" );
  iEvent.put(pf_mht_eta,  "pfmhteta" );
  iEvent.put(pf_mht_phi,  "pfmhtphi" );

  iEvent.put(calo_ht, "caloht");
  iEvent.put(calo_mht_pt, "calomhtpt");
  iEvent.put(calo_mht_eta, "calomhteta");
  iEvent.put(calo_mht_phi, "calomhtphi");

  iEvent.put(ngenlep,  "ngenlep" );
}

//define this as a plug-in
DEFINE_FWK_MODULE(babymaker);
