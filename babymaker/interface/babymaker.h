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
    edm::InputTag ElectronInputTag;
    edm::InputTag MuonInputTag;
    edm::InputTag pfJetsInputTag;
    edm::InputTag pfMetInputTag;
    edm::InputTag pfHTInputTag;
    edm::InputTag genJetsInputTag;
    
};


#endif
