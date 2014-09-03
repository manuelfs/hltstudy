import FWCore.ParameterSet.Config as cms

babymakerel = cms.EDProducer(
    "babymaker",
    MuonInputTag_       = cms.InputTag("unused"),
    #ElectronInputTag_   = cms.InputTag("hltPixelMatchElectronsActivity"), 
    ElectronInputTag_   = cms.InputTag("hltEgammaGsfElectrons"), 
    #pfJetsInputTag_     = cms.InputTag("hltAK4PFJetL1FastL2L3CorrectedNoPU"),
    pfJetsInputTag_     = cms.InputTag("hltAK4PFJetsCorrected"),
    pfMetInputTag_      = cms.InputTag("hltPFMETProducer"),
    #pfHTInputTag_       = cms.InputTag("hltEle15CaloIdTTrkIdTCaloIsoVLTrkIsoVLCleanedPFHTNoPU"),
    pfHTInputTag_       = cms.InputTag("hltPFHT"),
    genJetsInputTag_    = cms.InputTag("ak4GenJetsNoNu"),
)


