import FWCore.ParameterSet.Config as cms

babymaker = cms.EDProducer(
    "babymaker",
    hltMuonInputTag_       = cms.InputTag("unused"),
    #hltElectronInputTag_   = cms.InputTag("hltEgammaGsfElectrons"), 
    hltElectronInputTag_   = cms.InputTag("hltPixelMatchElectronsActivity"), 
    #hltPfJetsInputTag_     = cms.InputTag("hltAK4PFJetsCorrected"),
    hltPfJetsInputTag_     = cms.InputTag("hltAK4PFJetL1FastL2L3CorrectedNoPU"),
    hltPfMetInputTag_      = cms.InputTag("hltPFMETProducer"),
    #hltPfHTInputTag_       = cms.InputTag("hltPFHT"),
    hltPfHTInputTag_       = cms.InputTag("hltEle15CaloIdTTrkIdTCaloIsoVLTrkIsoVLCleanedPFHTNoPU"),
    hltGenJetsInputTag_    = cms.InputTag("ak4GenJetsNoNu"),
    recoMuonInputTag_      = cms.InputTag("muons"),
    recoElectronInputTag_  = cms.InputTag("gedGsfElectrons"), 
    recoPfJetsInputTag_    = cms.InputTag("ak4PFJets"),
    recoPfMetInputTag_     = cms.InputTag("pfMet"),
    recoGenJetsInputTag_   = cms.InputTag("ak5GenJets"),
)
