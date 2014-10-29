// IsoElectronProducer: Takes input from a filter and writes it as output to the event

#include "hltstudy/babymaker/interface/IsoElectronProducer.hpp"

#include <vector>
#include <memory>
#include <algorithm>
#include <limits>
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

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Math/interface/deltaR.h"

using  namespace std;
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
  track_map_token_(consumes<reco::RecoEcalCandidateIsolationMap>(track_map_tag_)),

  clustershape_tag_(iConfig.getParameter<edm::InputTag>("clustershapeTag")),
  clustershape_map_tag_(iConfig.getParameter<edm::InputTag>("clustershapeMapTag")),
  clustershape_token_(consumes<trigger::TriggerFilterObjectWithRefs>(clustershape_tag_)),
  clustershape_map_token_(consumes<reco::RecoEcalCandidateIsolationMap>(clustershape_map_tag_)),
  he_tag_(iConfig.getParameter<edm::InputTag>("heTag")),
  he_map_tag_(iConfig.getParameter<edm::InputTag>("heMapTag")),
  he_token_(consumes<trigger::TriggerFilterObjectWithRefs>(he_tag_)),
  he_map_token_(consumes<reco::RecoEcalCandidateIsolationMap>(he_map_tag_)),
  eminusp_tag_(iConfig.getParameter<edm::InputTag>("eminuspTag")),
  eminusp_map_tag_(iConfig.getParameter<edm::InputTag>("eminuspMapTag")),
  eminusp_token_(consumes<trigger::TriggerFilterObjectWithRefs>(eminusp_tag_)),
  eminusp_map_token_(consumes<reco::RecoEcalCandidateIsolationMap>(eminusp_map_tag_)),
  deta_tag_(iConfig.getParameter<edm::InputTag>("detaTag")),
  deta_map_tag_(iConfig.getParameter<edm::InputTag>("detaMapTag")),
  deta_token_(consumes<trigger::TriggerFilterObjectWithRefs>(deta_tag_)),
  deta_map_token_(consumes<reco::RecoEcalCandidateIsolationMap>(deta_map_tag_)),
  dphi_tag_(iConfig.getParameter<edm::InputTag>("dphiTag")),
  dphi_map_tag_(iConfig.getParameter<edm::InputTag>("dphiMapTag")),
  dphi_token_(consumes<trigger::TriggerFilterObjectWithRefs>(dphi_tag_)),
  dphi_map_token_(consumes<reco::RecoEcalCandidateIsolationMap>(dphi_map_tag_))
{

  produces<std::vector<float> >("elept").setBranchAlias("ele_pt");
  produces<std::vector<float> >("elephi").setBranchAlias("ele_phi");
  produces<std::vector<float> >("eleeta").setBranchAlias("ele_eta");
  produces<std::vector<float> >("eleecaliso").setBranchAlias("ele_ecal_iso"); 
  produces<std::vector<float> >("elehcaliso").setBranchAlias("ele_hcal_iso"); 
  produces<std::vector<float> >("eletrackiso").setBranchAlias("ele_track_iso"); 

  produces<std::vector<float> >("eleclustershape").setBranchAlias("ele_clustershape"); 

  produces<std::vector<float> >("elehe").setBranchAlias("ele_he"); 

  produces<std::vector<float> >("eleeminusp").setBranchAlias("ele_eminusp"); 

  produces<std::vector<float> >("eledeta").setBranchAlias("ele_deta"); 

  produces<std::vector<float> >("eledphi").setBranchAlias("ele_dphi"); 

}

IsoElectronProducer::~IsoElectronProducer(){
}

void IsoElectronProducer::beginJob(){
}

void IsoElectronProducer::endJob(){
}

void IsoElectronProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  auto ele_ecal_pt = make_auto(new std::vector<float>);
  auto ele_ecal_phi = make_auto(new std::vector<float>);
  auto ele_ecal_eta = make_auto(new std::vector<float>);
  auto ele_ecal_iso = make_auto(new std::vector<float>);
  auto ele_ecal_iso_temp = make_auto(new std::vector<float>);
  auto ele_hcal_pt = make_auto(new std::vector<float>);
  auto ele_hcal_phi = make_auto(new std::vector<float>);
  auto ele_hcal_eta = make_auto(new std::vector<float>);
  auto ele_hcal_iso = make_auto(new std::vector<float>);
  auto ele_hcal_iso_temp = make_auto(new std::vector<float>);
  auto ele_track_pt = make_auto(new std::vector<float>);
  auto ele_track_phi = make_auto(new std::vector<float>);
  auto ele_track_eta = make_auto(new std::vector<float>);
  auto ele_track_iso = make_auto(new std::vector<float>);

  auto ele_clustershape_pt = make_auto(new std::vector<float>);
  auto ele_clustershape_phi = make_auto(new std::vector<float>);
  auto ele_clustershape_eta = make_auto(new std::vector<float>);
  auto ele_clustershape_value = make_auto(new std::vector<float>);
  auto ele_clustershape_value_temp = make_auto(new std::vector<float>);
  auto ele_he_pt = make_auto(new std::vector<float>);
  auto ele_he_phi = make_auto(new std::vector<float>);
  auto ele_he_eta = make_auto(new std::vector<float>);
  auto ele_he_value = make_auto(new std::vector<float>);
  auto ele_he_value_temp = make_auto(new std::vector<float>);
  auto ele_eminusp_pt = make_auto(new std::vector<float>);
  auto ele_eminusp_phi = make_auto(new std::vector<float>);
  auto ele_eminusp_eta = make_auto(new std::vector<float>);
  auto ele_eminusp_value = make_auto(new std::vector<float>);
  auto ele_eminusp_value_temp = make_auto(new std::vector<float>);
  auto ele_deta_pt = make_auto(new std::vector<float>);
  auto ele_deta_phi = make_auto(new std::vector<float>);
  auto ele_deta_eta = make_auto(new std::vector<float>);
  auto ele_deta_value = make_auto(new std::vector<float>);
  auto ele_deta_value_temp = make_auto(new std::vector<float>);
  auto ele_dphi_pt = make_auto(new std::vector<float>);
  auto ele_dphi_phi = make_auto(new std::vector<float>);
  auto ele_dphi_eta = make_auto(new std::vector<float>);
  auto ele_dphi_value = make_auto(new std::vector<float>);
  auto ele_dphi_value_temp = make_auto(new std::vector<float>);
  
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
  track_iso_handle->getObjects(trigger::TriggerPhoton, track_ref);
  
  IsoMapLookup(ecal_ref, ecal_dep_map, ele_ecal_iso_temp, ele_ecal_pt, ele_ecal_phi, ele_ecal_eta);
  IsoMapLookup(hcal_ref, hcal_dep_map, ele_hcal_iso_temp, ele_hcal_pt, ele_hcal_phi, ele_hcal_eta);
  IsoMapLookup(track_ref, track_dep_map, ele_track_iso, ele_track_pt, ele_track_phi, ele_track_eta);

  GetSubset(ele_ecal_iso, ele_ecal_iso_temp, ele_ecal_pt, ele_track_pt);
  GetSubset(ele_hcal_iso, ele_hcal_iso_temp, ele_hcal_pt, ele_track_pt);

  edm::Handle<trigger::TriggerFilterObjectWithRefs> clustershape_value_handle;
  iEvent.getByToken(clustershape_token_, clustershape_value_handle);
  edm::Handle<reco::RecoEcalCandidateIsolationMap> clustershape_dep_map;
  iEvent.getByToken(clustershape_map_token_, clustershape_dep_map);
  std::vector<reco::RecoEcalCandidateRef> clustershape_ref;
  clustershape_value_handle->getObjects(trigger::TriggerCluster, clustershape_ref);
  IsoMapLookup(clustershape_ref, clustershape_dep_map, ele_clustershape_value_temp, 
	       ele_clustershape_pt, ele_clustershape_phi, ele_clustershape_eta, false, false);
  GetSubset(ele_clustershape_value, ele_clustershape_value_temp, ele_clustershape_pt, ele_track_pt);

  edm::Handle<trigger::TriggerFilterObjectWithRefs> he_value_handle;
  iEvent.getByToken(he_token_, he_value_handle);
  edm::Handle<reco::RecoEcalCandidateIsolationMap> he_dep_map;
  iEvent.getByToken(he_map_token_, he_dep_map);
  std::vector<reco::RecoEcalCandidateRef> he_ref;
  he_value_handle->getObjects(trigger::TriggerCluster, he_ref);
  IsoMapLookup(he_ref, he_dep_map, ele_he_value_temp, ele_he_pt, ele_he_phi, ele_he_eta, true, false);
  GetSubset(ele_he_value, ele_he_value_temp, ele_he_pt, ele_track_pt);

  edm::Handle<trigger::TriggerFilterObjectWithRefs> eminusp_value_handle;
  iEvent.getByToken(eminusp_token_, eminusp_value_handle);
  edm::Handle<reco::RecoEcalCandidateIsolationMap> eminusp_dep_map;
  iEvent.getByToken(eminusp_map_token_, eminusp_dep_map);
  std::vector<reco::RecoEcalCandidateRef> eminusp_ref;
  eminusp_value_handle->getObjects(trigger::TriggerCluster, eminusp_ref);
  IsoMapLookup(eminusp_ref, eminusp_dep_map, ele_eminusp_value_temp, ele_eminusp_pt, 
	       ele_eminusp_phi, ele_eminusp_eta, false);
  GetSubset(ele_eminusp_value, ele_eminusp_value_temp, ele_eminusp_pt, ele_track_pt);

  edm::Handle<trigger::TriggerFilterObjectWithRefs> deta_value_handle;
  iEvent.getByToken(deta_token_, deta_value_handle);
  edm::Handle<reco::RecoEcalCandidateIsolationMap> deta_dep_map;
  iEvent.getByToken(deta_map_token_, deta_dep_map);
  std::vector<reco::RecoEcalCandidateRef> deta_ref;
  deta_value_handle->getObjects(trigger::TriggerCluster, deta_ref);
  IsoMapLookup(deta_ref, deta_dep_map, ele_deta_value_temp, ele_deta_pt, ele_deta_phi, ele_deta_eta, false);
  GetSubset(ele_deta_value, ele_deta_value_temp, ele_deta_pt, ele_track_pt);

  edm::Handle<trigger::TriggerFilterObjectWithRefs> dphi_value_handle;
  iEvent.getByToken(dphi_token_, dphi_value_handle);
  edm::Handle<reco::RecoEcalCandidateIsolationMap> dphi_dep_map;
  iEvent.getByToken(dphi_map_token_, dphi_dep_map);
  std::vector<reco::RecoEcalCandidateRef> dphi_ref;
  dphi_value_handle->getObjects(trigger::TriggerCluster, dphi_ref);
  IsoMapLookup(dphi_ref, dphi_dep_map, ele_dphi_value_temp, ele_dphi_pt, ele_dphi_phi, ele_dphi_eta, false);
  GetSubset(ele_dphi_value, ele_dphi_value_temp, ele_dphi_pt, ele_track_pt);

  iEvent.put(ele_track_pt, "elept");
  iEvent.put(ele_track_phi, "elephi");
  iEvent.put(ele_track_eta, "eleeta");
  iEvent.put(ele_ecal_iso, "eleecaliso");
  iEvent.put(ele_hcal_iso, "elehcaliso");
  iEvent.put(ele_track_iso, "eletrackiso");

  iEvent.put(ele_clustershape_value, "eleclustershape");
  iEvent.put(ele_he_value, "elehe");
  iEvent.put(ele_eminusp_value, "eleeminusp");
  iEvent.put(ele_deta_value, "eledeta");
  iEvent.put(ele_dphi_value, "eledphi");



}

void GetSubset(std::auto_ptr<std::vector<float> >& iso,
	       std::auto_ptr<std::vector<float> >& temp_iso,
	       std::auto_ptr<std::vector<float> >& test_pt,
	       std::auto_ptr<std::vector<float> >& ref_pt){
  iso->resize(ref_pt->size());
  for(unsigned refi = 0; refi<ref_pt->size(); ++refi){
    bool found(false);
    for(unsigned testi = 0; testi<test_pt->size() && !found; ++testi){
      if(ref_pt->at(refi)==test_pt->at(testi)){
	iso->at(refi) = temp_iso->at(testi);
	found=true;
      }
    }
    if(!found){
      cms::Exception ex("MatchingFailure");
      ex.append("Could not match the tracker pt to any calo pt.");
      ex.addContext("Calling GetSubset");
      ex.addAdditionalInfo("IsoElectronProducer will now give up.");
      throw ex;
    }
  }
}

void IsoMapLookup(const std::vector<reco::RecoEcalCandidateRef>& ref,
		  const edm::Handle<reco::RecoEcalCandidateIsolationMap>& map,
		  std::auto_ptr<std::vector<float> >& iso,
		  std::auto_ptr<std::vector<float> >& pt,
		  std::auto_ptr<std::vector<float> >& phi,
		  std::auto_ptr<std::vector<float> >& eta, bool overE, bool useEt){
  if(iso.get()) iso->clear();
  if(pt.get()) pt->clear();
  if(phi.get()) phi->clear();
  if(eta.get()) eta->clear();
  
  for(const auto& ele: ref){
    const float etaSC = ele->eta();
    float energy = ele->superCluster()->energy();
    if(useEt) energy *= sin(2.0*atan(exp(-etaSC)));
    float vali = map->find(ele)->val;
    if(overE) vali /= energy;
      
    if (iso.get()) iso->push_back(vali);
    if (pt.get()) pt->push_back(ele->pt());
    if (phi.get()) phi->push_back(ele->phi());
    if (eta.get()) eta->push_back(ele->eta());
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoElectronProducer);
