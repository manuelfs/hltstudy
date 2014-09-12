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
};

void IsoMapLookup(const std::vector<reco::RecoEcalCandidateRef>& ref,
		  const edm::Handle<reco::RecoEcalCandidateIsolationMap>& map,
		  std::auto_ptr<std::vector<float> >& iso,
		  std::auto_ptr<std::vector<float> >& pt,
		  std::auto_ptr<std::vector<float> >& phi,
		  std::auto_ptr<std::vector<float> >& eta);

void GetSubset(std::auto_ptr<std::vector<float> >& iso,
	       std::auto_ptr<std::vector<float> >& temp_iso,
	       std::auto_ptr<std::vector<float> >& test_pt,
	       std::auto_ptr<std::vector<float> >& ref_pt);

#endif
