// IsoElectronProducer: Takes input from a filter and writes it as output to the event

#include "hltstudy/babymaker/interface/IsoElectronProducer.hpp"

#include <vector>
#include <iostream>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

namespace{
  template<typename T>
    std::auto_ptr<T> make_auto(T * const p=nullptr){
    return std::auto_ptr<T>{p};
  }
}

IsoElectronProducer::IsoElectronProducer(const edm::ParameterSet& iConfig):
  ecal_tag_(iConfig.getParameter<edm::InputTag>("ecalTag")),
  hcal_tag_(iConfig.getParameter<edm::InputTag>("hcalTag")),
  track_tag_(iConfig.getParameter<edm::InputTag>("trackTag")),
  ecal_map_tag_(iConfig.getParameter<edm::InputTag>("ecalMapTag")),
  hcal_map_tag_(iConfig.getParameter<edm::InputTag>("hcalMapTag")),
  track_map_tag_(iConfig.getParameter<edm::InputTag>("trackMapTag")),
  ecal_token_(consumes<trigger::TriggerFilterObjectWithRefs>(ecal_tag_)),
  hcal_token_(consumes<trigger::TriggerFilterObjectWithRefs>(hcal_tag_)),
  track_token_(consumes<trigger::TriggerFilterObjectWithRefs>(track_tag_)),
  ecal_map_token_(consumes<reco::RecoEcalCandidateIsolationMap>(ecal_map_tag_)),
  hcal_map_token_(consumes<reco::RecoEcalCandidateIsolationMap>(hcal_map_tag_)),
  track_map_token_(consumes<reco::RecoEcalCandidateIsolationMap>(track_map_tag_)){
  produces<std::vector<float> >("elept").setBranchAlias("ele_pt");
  produces<std::vector<float> >("elephi").setBranchAlias("ele_phi");
  produces<std::vector<float> >("eleeta").setBranchAlias("ele_eta");
  produces<std::vector<float> >("eleecaliso").setBranchAlias("ele_ecal_iso"); 
  produces<std::vector<float> >("elehcaliso").setBranchAlias("ele_hcal_iso"); 
  produces<std::vector<float> >("eletrackiso").setBranchAlias("ele_track_iso"); 
}

IsoElectronProducer::~IsoElectronProducer(){
}

void IsoElectronProducer::beginJob(){
}

void IsoElectronProducer::endJob(){
}

void IsoElectronProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  auto ele_pt = make_auto(new std::vector<float>);
  auto ele_phi = make_auto(new std::vector<float>);
  auto ele_eta = make_auto(new std::vector<float>);
  auto ele_ecal_iso = make_auto(new std::vector<float>);
  auto ele_hcal_iso = make_auto(new std::vector<float>);
  auto ele_track_iso = make_auto(new std::vector<float>);
  
  edm::Handle<trigger::TriggerFilterObjectWithRefs> ecal_iso_handle;
  edm::Handle<trigger::TriggerFilterObjectWithRefs> hcal_iso_handle;
  edm::Handle<trigger::TriggerFilterObjectWithRefs> track_iso_handle;
  iEvent.getByToken(ecal_token_, ecal_iso_handle);
  iEvent.getByToken(hcal_token_, hcal_iso_handle);
  iEvent.getByToken(track_token_, track_iso_handle);
  
  edm::Handle<reco::RecoEcalCandidateIsolationMap> ecal_dep_map;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> hcal_dep_map;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> track_dep_map;
  iEvent.getByToken(ecal_map_token_, ecal_dep_map);
  iEvent.getByToken(hcal_map_token_, hcal_dep_map);
  iEvent.getByToken(track_map_token_, track_dep_map);
  
  std::vector<reco::RecoEcalCandidateRef> ecal_ref;
  ecal_iso_handle->getObjects(trigger::TriggerCluster, ecal_ref);
  std::vector<reco::RecoEcalCandidateRef> hcal_ref;
  hcal_iso_handle->getObjects(trigger::TriggerCluster, hcal_ref);
  std::vector<reco::RecoEcalCandidateRef> track_ref;
  track_iso_handle->getObjects(trigger::TriggerCluster, track_ref);
  
  for(std::vector<reco::RecoEcalCandidateRef>::const_iterator it = ecal_ref.begin(); it != ecal_ref.end(); ++it){
    reco::RecoEcalCandidateIsolationMap::const_iterator ecal_mapi = ecal_dep_map->find(*it);
    float vali = ecal_mapi->val;
    float energy = (*it)->superCluster()->energy();
    float etaSC = (*it)->eta();
    energy *= sin(2.0 * atan(exp(-etaSC)));
    ele_ecal_iso->push_back(vali/energy);
    
    ele_pt->push_back((*it)->pt());
    ele_phi->push_back((*it)->phi());
    ele_eta->push_back((*it)->eta());
  }

  for(std::vector<reco::RecoEcalCandidateRef>::const_iterator it = hcal_ref.begin(); it != hcal_ref.end(); ++it){
    reco::RecoEcalCandidateIsolationMap::const_iterator hcal_mapi = hcal_dep_map->find(*it);
    float vali = hcal_mapi->val;
    float energy = (*it)->superCluster()->energy();
    float etaSC = (*it)->eta();
    energy *= sin(2.0 * atan(exp(-etaSC)));
    ele_hcal_iso->push_back(vali/energy);
  }

  for(std::vector<reco::RecoEcalCandidateRef>::const_iterator it = track_ref.begin(); it != track_ref.end(); ++it){
    reco::RecoEcalCandidateIsolationMap::const_iterator track_mapi = track_dep_map->find(*it);
    float vali = track_mapi->val;
    float energy = (*it)->superCluster()->energy();
    float etaSC = (*it)->eta();
    energy *= sin(2.0 * atan(exp(-etaSC)));
    ele_track_iso->push_back(vali/energy);
  }
  
  if(ele_ecal_iso->size()<ele_pt->size()) ele_ecal_iso->resize(ele_pt->size(), -1.0);
  if(ele_hcal_iso->size()<ele_pt->size()) ele_hcal_iso->resize(ele_pt->size(), -1.0);
  if(ele_track_iso->size()<ele_pt->size()) ele_track_iso->resize(ele_pt->size(), -1.0);
  
  iEvent.put(ele_pt, "elept");
  iEvent.put(ele_phi, "elephi");
  iEvent.put(ele_eta, "eleeta");
  iEvent.put(ele_ecal_iso, "eleecaliso");
  iEvent.put(ele_hcal_iso, "elehcaliso");
  iEvent.put(ele_track_iso, "eletrackiso");
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoElectronProducer);
