import FWCore.ParameterSet.Config as cms

babymaker = cms.EDProducer(
    "babymaker",
    hltGenJetsInputTag_    = cms.InputTag("ak4GenJetsNoNu"),
    hltMuonInputTag_       = cms.InputTag("unused"),
    hltElectronInputTag_   = cms.InputTag("unused"), 
    hltPfJetsInputTag_     = cms.InputTag("unused"),
    hltPfMetInputTag_      = cms.InputTag("unused"),
    hltPfHTInputTag_       = cms.InputTag("unused"),
    recoMuonInputTag_      = cms.InputTag("unused"),
    recoElectronInputTag_  = cms.InputTag("unused"), 
    recoPfJetsInputTag_    = cms.InputTag("unused"),
    recoPfMetInputTag_     = cms.InputTag("unused"),
    recoGenJetsInputTag_   = cms.InputTag("unused"),
)
