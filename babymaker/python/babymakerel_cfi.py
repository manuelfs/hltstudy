import FWCore.ParameterSet.Config as cms

babymakerel = cms.EDProducer(
    "babymaker",
    MuonInputTag_       = cms.InputTag("unused"),
    ElectronInputTag_   = cms.InputTag("hltPixelMatchElectronsActivity"), 
    pfJetsInputTag_     = cms.InputTag("hltAK4PFJetL1FastL2L3CorrectedNoPU"),
    pfMetInputTag_      = cms.InputTag("hltPFMETProducer"),
    pfHTInputTag_       = cms.InputTag("hltEle15CaloIdTTrkIdTCaloIsoVLTrkIsoVLCleanedPFHTNoPU"),
    genJetsInputTag_    = cms.InputTag("ak4GenJetsNoNu"),
)


