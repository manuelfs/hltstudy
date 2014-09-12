import FWCore.ParameterSet.Config as cms

babymaker = cms.EDProducer(
    "babymaker",
    hltMuonInputString_     = cms.string("unused"),
    hltElectronInputString_ = cms.string("unused"),
    hltPfJetsInputTag_      = cms.InputTag("unused"),
    hltPfMetInputTag_       = cms.InputTag("unused"),
    hltPfHTInputTag_        = cms.InputTag("unused"),
    hltGenJetsInputTag_     = cms.InputTag("ak4GenJetsNoNu"),
    recoMuonInputTag_       = cms.InputTag("unused"),
    recoElectronInputTag_   = cms.InputTag("unused"), 
    recoPfJetsInputTag_     = cms.InputTag("unused"),
    recoPfMetInputTag_      = cms.InputTag("unused"),
    recoGenJetsInputTag_    = cms.InputTag("unused"),
    m_Jets_                 = cms.InputTag("unused"),
    m_JetTags_              = cms.InputTag("unused")
)
