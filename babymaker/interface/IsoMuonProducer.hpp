#ifndef __H_MU_PROD__
#define __H_MU_PROD__

#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

class IsoMuonProducer : public edm::EDProducer {
public:
  explicit IsoMuonProducer(const edm::ParameterSet&);
  ~IsoMuonProducer();

private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::InputTag iso_tag_;
  edm::InputTag cand_tag_;

  edm::EDGetTokenT<edm::ValueMap<double> > iso_token_;
  edm::EDGetTokenT<std::vector<reco::RecoChargedCandidate> > cand_token_;
};

#endif
