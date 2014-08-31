// babymaker: Makes flat ntuples with HT, MET, MC, jets, and leptons

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "hltstudy/babymaker/interface/babymaker.h"

#include "TMath.h"

babymaker::babymaker(const edm::ParameterSet& iConfig) {

  produces<float> ("wl1ht200").setBranchAlias("wl1ht200");
  produces<float> ("pfht").setBranchAlias("pf_ht");
  produces<float> ("genht").setBranchAlias("gen_ht");

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

  produces<float> ("metpt").setBranchAlias("met_pt");
  produces<float> ("meteta").setBranchAlias("met_eta");
  produces<float> ("metphi").setBranchAlias("met_phi");

  produces<float> ("pfmhtpt").setBranchAlias("pf_mht_pt");
  produces<float> ("pfmhteta").setBranchAlias("pf_mht_eta");
  produces<float> ("pfmhtphi").setBranchAlias("pf_mht_phi");

  produces<unsigned int> ("ngenlep").setBranchAlias("ngenlep");

  pfJetsInputTag = iConfig.getParameter<edm::InputTag>("pfJetsInputTag_");
  ElectronInputTag = iConfig.getParameter<edm::InputTag>("ElectronInputTag_");
  MuonInputTag = iConfig.getParameter<edm::InputTag>("MuonInputTag_");
  pfMetInputTag = iConfig.getParameter<edm::InputTag>("pfMetInputTag_");
  pfHTInputTag = iConfig.getParameter<edm::InputTag>("pfHTInputTag_");
  genJetsInputTag = iConfig.getParameter<edm::InputTag>("genJetsInputTag_");
}


babymaker::~babymaker() {}

void  babymaker::beginJob() {
}

void babymaker::endJob() {
}


// ------------ method called to produce the data  ------------
void babymaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
  std::auto_ptr<std::vector<float> > els_pt   (new std::vector<float>);
  std::auto_ptr<std::vector<float> > els_eta  (new std::vector<float>);
  std::auto_ptr<std::vector<float> > els_phi  (new std::vector<float>);
  
  std::auto_ptr<std::vector<float> > mus_pt   (new std::vector<float>);
  std::auto_ptr<std::vector<float> > mus_eta  (new std::vector<float>);
  std::auto_ptr<std::vector<float> > mus_phi  (new std::vector<float>);
    
  std::auto_ptr<std::vector<float> > pfjets_pt   (new std::vector<float>);
  std::auto_ptr<std::vector<float> > pfjets_eta  (new std::vector<float>);
  std::auto_ptr<std::vector<float> > pfjets_phi  (new std::vector<float>);

  std::auto_ptr<std::vector<float> > genjets_pt   (new std::vector<float>);
  std::auto_ptr<std::vector<float> > genjets_eta  (new std::vector<float>);
  std::auto_ptr<std::vector<float> > genjets_phi  (new std::vector<float>);

  std::auto_ptr<std::vector<float> > genels_pt   (new std::vector<float>);
  std::auto_ptr<std::vector<float> > genels_eta  (new std::vector<float>);
  std::auto_ptr<std::vector<float> > genels_phi  (new std::vector<float>);

  std::auto_ptr<std::vector<float> > genmus_pt   (new std::vector<float>);
  std::auto_ptr<std::vector<float> > genmus_eta  (new std::vector<float>);
  std::auto_ptr<std::vector<float> > genmus_phi  (new std::vector<float>);

  std::auto_ptr<float> met_pt   (new float);
  std::auto_ptr<float> met_eta  (new float);
  std::auto_ptr<float> met_phi  (new float);

  std::auto_ptr<float> pf_ht   (new float);
  std::auto_ptr<float> gen_ht   (new float);
  std::auto_ptr<float> wl1ht200   (new float);

  std::auto_ptr<float> pf_mht_pt   (new float);
  std::auto_ptr<float> pf_mht_eta   (new float);
  std::auto_ptr<float> pf_mht_phi   (new float);

  std::auto_ptr<float> pseudoMT2_hlt  (new float);
  std::auto_ptr<float> pseudoMT2_snt  (new float);
  std::auto_ptr<float> mt2_hlt  (new float);
  std::auto_ptr<float> mt2_snt  (new float);

  std::auto_ptr<unsigned int> ngenlep (new unsigned int);

  edm::Handle<edm::View<reco::PFJet> > jet_h;
  iEvent.getByLabel(pfJetsInputTag, jet_h);

  edm::Handle<edm::View<reco::MET> > pfht_h;
  iEvent.getByLabel(pfHTInputTag, pfht_h);

  edm::Handle<edm::View<reco::MET> > met_h;
  iEvent.getByLabel(pfMetInputTag, met_h);

  edm::Handle<edm::View<reco::GenJet> > genjet_h;
  iEvent.getByLabel(genJetsInputTag, genjet_h);

  edm::Handle<reco::GenParticleCollection> genparts;
  iEvent.getByLabel("genParticles", genparts);


  *ngenlep = 0;
  for (reco::GenParticleCollection::const_iterator genp = genparts->begin(); genp != genparts->end(); ++ genp ) {
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
  if(ElectronInputTag.label() !="unused") {
    iEvent.getByLabel(ElectronInputTag, els_h);
    for(edm::View<reco::Electron>::const_iterator els_it = els_h->begin(); els_it != els_h->end(); els_it++){
      if(abs(els_it->pdgId()) != 11) continue;
      els_pt  ->push_back(els_it->pt());
      els_eta ->push_back(els_it->eta());
      els_phi ->push_back(els_it->phi());
    } 
  }

  edm::Handle<edm::View<reco::Muon> > mus_h;
  if(MuonInputTag.label() !="unused") {
    iEvent.getByLabel(MuonInputTag, mus_h);
    for(edm::View<reco::Muon>::const_iterator mus_it = mus_h->begin(); mus_it != mus_h->end(); mus_it++){
      mus_pt  ->push_back(mus_it->pt());
      mus_eta ->push_back(mus_it->eta());
      mus_phi ->push_back(mus_it->phi());
    } 
  }

  for(edm::View<reco::PFJet>::const_iterator jet_it = jet_h->begin(); jet_it != jet_h->end(); jet_it++){

    if(jet_it->pt() < 35.0) continue;
    pfjets_pt  ->push_back(jet_it->pt());
    pfjets_eta ->push_back(jet_it->eta());
    pfjets_phi ->push_back(jet_it->phi());
  } 
  *gen_ht = 0;
  for(edm::View<reco::GenJet>::const_iterator genjet_it = genjet_h->begin(); genjet_it != genjet_h->end(); genjet_it++){

    if(genjet_it->pt() < 35.0) continue;

    genjets_pt  ->push_back(genjet_it->pt());
    genjets_eta ->push_back(genjet_it->eta());
    genjets_phi ->push_back(genjet_it->phi());

    if(genjet_it->pt()>40 && genjet_it->eta() < 3) *gen_ht += genjet_it->pt();

  } 
  // Turn on curve for L1 HTT200 (calculated by Dominick)
  *wl1ht200 = (0.5*TMath::Erf((1.35121e-02)*(*gen_ht-(3.02695e+02)))+0.5);

  *met_pt   = (met_h->front()).pt();
  *met_eta  = (met_h->front()).eta();
  *met_phi  = (met_h->front()).phi();

  *pf_ht    = (pfht_h->front()).sumEt();

  *pf_mht_pt    = (pfht_h->front()).pt();
  *pf_mht_eta   = (pfht_h->front()).eta();
  *pf_mht_phi   = (pfht_h->front()).phi();

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

  iEvent.put(met_pt,   "metpt" );
  iEvent.put(met_eta,  "meteta" );
  iEvent.put(met_phi,  "metphi" );

  iEvent.put(pf_ht,   "pfht" );
  iEvent.put(gen_ht,   "genht" );
  iEvent.put(wl1ht200,   "wl1ht200" );

  iEvent.put(pf_mht_pt,   "pfmhtpt" );
  iEvent.put(pf_mht_eta,  "pfmhteta" );
  iEvent.put(pf_mht_phi,  "pfmhtphi" );

  iEvent.put(ngenlep,  "ngenlep" );
}

//define this as a plug-in
DEFINE_FWK_MODULE(babymaker);
