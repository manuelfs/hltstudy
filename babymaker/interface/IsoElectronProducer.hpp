#ifndef __H_MU_PROD__
#define __H_MU_PROD__

#include <vector>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

class IsoElectronProducer : public edm::EDProducer {
public:
  explicit IsoElectronProducer(const edm::ParameterSet&);
  ~IsoElectronProducer();

private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::InputTag ecal_tag_;
  edm::InputTag hcal_tag_;
  edm::InputTag track_tag_;
  edm::InputTag ecal_map_tag_;
  edm::InputTag hcal_map_tag_;
  edm::InputTag track_map_tag_;

  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> ecal_token_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> hcal_token_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> track_token_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> ecal_map_token_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> hcal_map_token_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> track_map_token_;

  edm::InputTag clustershape_tag_;
  edm::InputTag clustershape_map_tag_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> clustershape_token_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> clustershape_map_token_;

  edm::InputTag he_tag_;
  edm::InputTag he_map_tag_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> he_token_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> he_map_token_;

  edm::InputTag eminusp_tag_;
  edm::InputTag eminusp_map_tag_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> eminusp_token_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> eminusp_map_token_;

  edm::InputTag deta_tag_;
  edm::InputTag deta_map_tag_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> deta_token_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> deta_map_token_;

  edm::InputTag dphi_tag_;
  edm::InputTag dphi_map_tag_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> dphi_token_;
  edm::EDGetTokenT<reco::RecoEcalCandidateIsolationMap> dphi_map_token_;


};

void IsoMapLookup(const std::vector<reco::RecoEcalCandidateRef>& ref,
		  const edm::Handle<reco::RecoEcalCandidateIsolationMap>& map,
		  std::auto_ptr<std::vector<float> >& iso,
		  std::auto_ptr<std::vector<float> >& pt,
		  std::auto_ptr<std::vector<float> >& phi,
		  std::auto_ptr<std::vector<float> >& eta, 
		  bool overE=true, bool useEt=true);

void GetSubset(std::auto_ptr<std::vector<float> >& iso,
	       std::auto_ptr<std::vector<float> >& temp_iso,
	       std::auto_ptr<std::vector<float> >& test_pt,
	       std::auto_ptr<std::vector<float> >& ref_pt);

#endif
