// IsoMuonProducer: Takes input from a filter and writes it as output to the event

#include "hltstudy/babymaker/interface/IsoMuonProducer.hpp"

#include <vector>
#include <iostream>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

namespace{
  template<typename T>
    std::auto_ptr<T> make_auto(T * const p=nullptr){
    return std::auto_ptr<T>{p};
  }
}

IsoMuonProducer::IsoMuonProducer(const edm::ParameterSet& iConfig):
  iso_tag_(iConfig.getParameter<edm::InputTag>("isoTag")),
  cand_tag_(iConfig.getParameter<edm::InputTag>("candTag")),
  iso_token_(consumes<edm::ValueMap<double> >(iso_tag_)),
  cand_token_(consumes<std::vector<reco::RecoChargedCandidate> >(cand_tag_)){
  produces<std::vector<float> >("mupt").setBranchAlias("mu_pt");
  produces<std::vector<float> >("muphi").setBranchAlias("mu_phi");
  produces<std::vector<float> >("mueta").setBranchAlias("mu_eta");
  produces<std::vector<float> >("muiso").setBranchAlias("mu_iso"); 
}

IsoMuonProducer::~IsoMuonProducer(){
}

void IsoMuonProducer::beginJob(){
}

void IsoMuonProducer::endJob(){
}

void IsoMuonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  auto mu_pt = make_auto(new std::vector<float>);
  auto mu_phi = make_auto(new std::vector<float>);
  auto mu_eta = make_auto(new std::vector<float>);
  auto mu_iso = make_auto(new std::vector<float>);
  
  edm::Handle<edm::ValueMap<double> > iso_handle;
  edm::Handle<std::vector<reco::RecoChargedCandidate> > muon_handle;
  iEvent.getByToken(iso_token_, iso_handle);
  iEvent.getByToken(cand_token_, muon_handle);
  //iEvent.getByLabel("hltL3MuonCombRelIsolationsIterTrkRegIter02","combinedRelativeIsoDeposits", iso_handle);
  //iEvent.getByLabel("hltL3MuonCandidates", muon_handle);
  
  for(auto it = iso_handle->begin(); it != iso_handle->end(); ++it){
    for(auto ite = it.begin(); ite != it.end(); ++ite){
      mu_iso->push_back(*ite);
    }
  }
  
  for(auto it = muon_handle->begin(); it != muon_handle->end(); ++it){
    mu_pt->push_back(it->pt());
    mu_phi->push_back(it->phi());
    mu_eta->push_back(it->eta());
  }

  iEvent.put(mu_pt, "mupt");
  iEvent.put(mu_phi, "muphi");
  iEvent.put(mu_eta, "mueta");
  iEvent.put(mu_iso, "muiso");
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoMuonProducer);
