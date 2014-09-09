#ifndef babymaker_H
#define babymaker_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class babymaker : public edm::EDProducer {
 public:
  explicit babymaker (const edm::ParameterSet&);
  ~babymaker();
  
 private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  edm::InputTag hltElectronInputTag;
  edm::InputTag hltMuonInputTag;
  edm::InputTag hltPfJetsInputTag;
  edm::InputTag hltPfMetInputTag;
  edm::InputTag hltPfHTInputTag;
  edm::InputTag hltGenJetsInputTag;
  
  edm::InputTag recoElectronInputTag;
  edm::InputTag recoMuonInputTag;
  edm::InputTag recoPfJetsInputTag;
  edm::InputTag recoPfMetInputTag;
  edm::InputTag recoGenJetsInputTag;
  edm::InputTag m_Jets;   
  edm::InputTag m_JetTags; 
  
  edm::EDGetTokenT<std::vector<reco::CaloJet> > m_JetsToken; 
  edm::EDGetTokenT<reco::JetTagCollection> m_JetTagsToken;
};


#endif
