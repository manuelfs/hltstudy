import FWCore.ParameterSet.Config as cms

babymakermu = cms.EDProducer(
    "babymaker",
    MuonInputTag_       = cms.InputTag("hltMuons"),
    ElectronInputTag_   = cms.InputTag("unused"), 
    pfJetsInputTag_     = cms.InputTag("hltAK4PFJetL1FastL2L3CorrectedNoPU"),
    pfMetInputTag_      = cms.InputTag("hltPFMETProducer"),
    pfHTInputTag_       = cms.InputTag("hltPFHTNoPU"),
    genJetsInputTag_    = cms.InputTag("ak4GenJetsNoNu"),
)
