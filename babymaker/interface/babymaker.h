#ifndef babymaker_H
#define babymaker_H

// system include files
#include <string>
#include <memory>
#include <vector>

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
  std::string hltElectronInputString;
  edm::InputTag hltElectronPtInputTag;
  edm::InputTag hltElectronPhiInputTag;
  edm::InputTag hltElectronEtaInputTag;
  edm::InputTag hltElectronEcalIsoInputTag;
  edm::InputTag hltElectronHcalIsoInputTag;
  edm::InputTag hltElectronTrackIsoInputTag;
  std::string hltMuonInputString;
  edm::InputTag hltMuonPtInputTag;
  edm::InputTag hltMuonPhiInputTag;
  edm::InputTag hltMuonEtaInputTag;
  edm::InputTag hltMuonIsoInputTag;
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
