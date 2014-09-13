import FWCore.ParameterSet.Config as cms

babymaker = cms.EDProducer(
    "babymaker",
    hltMuonInputString_     = cms.string("IsoMuonProducer"),
    hltElectronInputString_ = cms.string("unused"),
    hltPfJetsInputTag_      = cms.InputTag("hltAK4PFJetsCorrected"),
    hltPfMetInputTag_       = cms.InputTag("hltPFMETProducer"),
    hltPfHTInputTag_        = cms.InputTag("hltPFHT"),
    hltCaloMetInputTag_     = cms.InputTag("hltHtMht"),
    hltCaloHTInputTag_      = cms.InputTag("hltHtMht"),
    hltGenJetsInputTag_     = cms.InputTag("ak4GenJetsNoNu"),
    recoMuonInputTag_       = cms.InputTag("unused"),
    recoElectronInputTag_   = cms.InputTag("unused"), 
    recoPfJetsInputTag_     = cms.InputTag("unused"),
    recoPfMetInputTag_      = cms.InputTag("unused"),
    recoGenJetsInputTag_    = cms.InputTag("unused"),
    m_Jets_                 = cms.InputTag("hltSelector4CentralJetsL1FastJet"), # btagging 
    m_JetTags_              = cms.InputTag("hltL3CombinedSecondaryVertexBJetTags"), # btagging 
)
