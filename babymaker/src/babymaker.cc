// babymaker: Makes flat ntuples with HT, MET, MC, jets, and leptons

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
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
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

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
  produces<float> ("genht").setBranchAlias("gen_ht");

  produces<float> ("genmet").setBranchAlias("gen_met");
  produces<float> ("genmetcalo").setBranchAlias("gen_metcalo");
  produces<float> ("genmetcalononprompt").setBranchAlias("gen_metcalononprompt");

  produces<std::vector<float> > ("elspt").setBranchAlias("els_pt");
  produces<std::vector<float> > ("elseta").setBranchAlias("els_eta");
  produces<std::vector<float> > ("elsphi").setBranchAlias("els_phi");
    
  produces<std::vector<float> > ("muspt").setBranchAlias("mus_pt");
  produces<std::vector<float> > ("museta").setBranchAlias("mus_eta");
  produces<std::vector<float> > ("musphi").setBranchAlias("mus_phi");
    
  produces<std::vector<float> > ("pfjetspt").setBranchAlias("pfjets_pt");
  produces<std::vector<float> > ("pfjetseta").setBranchAlias("pfjets_eta");
  produces<std::vector<float> > ("pfjetsphi").setBranchAlias("pfjets_phi");

  produces<std::vector<float> > ("genjetspt").setBranchAlias("genjets_pt");
  produces<std::vector<float> > ("genjetseta").setBranchAlias("genjets_eta");
  produces<std::vector<float> > ("genjetsphi").setBranchAlias("genjets_phi");
    
  produces<std::vector<float> > ("genelspt").setBranchAlias("genels_pt");
  produces<std::vector<float> > ("genelseta").setBranchAlias("genels_eta");
  produces<std::vector<float> > ("genelsphi").setBranchAlias("genels_phi");
    
  produces<std::vector<float> > ("genmuspt").setBranchAlias("genmus_pt");
  produces<std::vector<float> > ("genmuseta").setBranchAlias("genmus_eta");
  produces<std::vector<float> > ("genmusphi").setBranchAlias("genmus_phi");

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

  produces<float> ("metpt").setBranchAlias("met_pt");
  produces<float> ("meteta").setBranchAlias("met_eta");
  produces<float> ("metphi").setBranchAlias("met_phi");

  produces<float> ("pfmhtpt").setBranchAlias("pf_mht_pt");
  produces<float> ("pfmhteta").setBranchAlias("pf_mht_eta");
  produces<float> ("pfmhtphi").setBranchAlias("pf_mht_phi");

  produces<unsigned int> ("ngenlep").setBranchAlias("ngenlep");

  hltPfJetsInputTag = iConfig.getParameter<edm::InputTag>("hltPfJetsInputTag_");
  hltElectronInputTag = iConfig.getParameter<edm::InputTag>("hltElectronInputTag_");
  hltMuonInputTag = iConfig.getParameter<edm::InputTag>("hltMuonInputTag_");
  hltPfMetInputTag = iConfig.getParameter<edm::InputTag>("hltPfMetInputTag_");
  hltPfHTInputTag = iConfig.getParameter<edm::InputTag>("hltPfHTInputTag_");
  hltGenJetsInputTag = iConfig.getParameter<edm::InputTag>("hltGenJetsInputTag_");

  recoPfJetsInputTag = iConfig.getParameter<edm::InputTag>("recoPfJetsInputTag_");
  recoElectronInputTag = iConfig.getParameter<edm::InputTag>("recoElectronInputTag_");
  recoMuonInputTag = iConfig.getParameter<edm::InputTag>("recoMuonInputTag_");
  recoPfMetInputTag = iConfig.getParameter<edm::InputTag>("recoPfMetInputTag_");
  recoGenJetsInputTag = iConfig.getParameter<edm::InputTag>("recoGenJetsInputTag_");
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
  
  auto mus_pt  = make_auto(new std::vector<float>);
  auto mus_eta = make_auto(new std::vector<float>);
  auto mus_phi = make_auto(new std::vector<float>);
    
  auto pfjets_pt  = make_auto(new std::vector<float>);
  auto pfjets_eta = make_auto(new std::vector<float>);
  auto pfjets_phi = make_auto(new std::vector<float>);

  auto genjets_pt  = make_auto(new std::vector<float>);
  auto genjets_eta = make_auto(new std::vector<float>);
  auto genjets_phi = make_auto(new std::vector<float>);

  auto genels_pt  = make_auto(new std::vector<float>);
  auto genels_eta = make_auto(new std::vector<float>);
  auto genels_phi = make_auto(new std::vector<float>);

  auto genmus_pt  = make_auto(new std::vector<float>);
  auto genmus_eta = make_auto(new std::vector<float>);
  auto genmus_phi = make_auto(new std::vector<float>);

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
        
    if (genp->status() != 3 && genp->status() != 23) continue;  

    int id = genp->pdgId();
        
    if( TMath::Abs(id) != 11 && TMath::Abs(id) != 13 ) continue; 
        
    if( TMath::Abs(id) == 11 ) {
      genels_pt  ->push_back(genp->p4().pt());
      genels_eta ->push_back(genp->p4().eta());
      genels_phi ->push_back(genp->p4().phi());
    }
    if( TMath::Abs(id) == 13 ) {
      genmus_pt  ->push_back(genp->p4().pt());
      genmus_eta ->push_back(genp->p4().eta());
      genmus_phi ->push_back(genp->p4().phi());
    }
  }
    

  edm::Handle<edm::View<reco::Electron> > els_h;
  if(hltElectronInputTag.label() !="unused") {
    iEvent.getByLabel(hltElectronInputTag, els_h);
    for(auto els_it = els_h->begin(); els_it != els_h->end(); ++els_it){
      if(abs(els_it->pdgId()) != 11) continue;
      els_pt  ->push_back(els_it->pt());
      els_eta ->push_back(els_it->eta());
      els_phi ->push_back(els_it->phi());
    } 
  }

  edm::Handle<edm::View<reco::Muon> > mus_h;
  if(hltMuonInputTag.label() !="unused") {
    iEvent.getByLabel(hltMuonInputTag, mus_h);
    for(auto mus_it = mus_h->begin(); mus_it != mus_h->end(); ++mus_it){
      mus_pt  ->push_back(mus_it->pt());
      mus_eta ->push_back(mus_it->eta());
      mus_phi ->push_back(mus_it->phi());
    } 
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
  
  // Turn on curve for L1 HTT200 (calculated by Dominick)
  *wl1ht200 = (0.5*TMath::Erf((1.35121e-02)*(*gen_ht-(3.02695e+02)))+0.5);



  iEvent.put(els_pt,   "elspt" );
  iEvent.put(els_eta,  "elseta" );
  iEvent.put(els_phi,  "elsphi" );
    
  iEvent.put(mus_pt,   "muspt" );
  iEvent.put(mus_eta,  "museta" );
  iEvent.put(mus_phi,  "musphi" );

  iEvent.put(pfjets_pt,   "pfjetspt" );
  iEvent.put(pfjets_eta,  "pfjetseta" );
  iEvent.put(pfjets_phi,  "pfjetsphi" );

  iEvent.put(genjets_pt,   "genjetspt" );
  iEvent.put(genjets_eta,  "genjetseta" );
  iEvent.put(genjets_phi,  "genjetsphi" );

  iEvent.put(genels_pt,   "genelspt" );
  iEvent.put(genels_eta,  "genelseta" );
  iEvent.put(genels_phi,  "genelsphi" );

  iEvent.put(genmus_pt,   "genmuspt" );
  iEvent.put(genmus_eta,  "genmuseta" );
  iEvent.put(genmus_phi,  "genmusphi" );

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

  iEvent.put(ngenlep,  "ngenlep" );
}

//define this as a plug-in
DEFINE_FWK_MODULE(babymaker);
