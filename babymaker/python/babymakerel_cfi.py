import FWCore.ParameterSet.Config as cms

babymaker = cms.EDProducer(
    "babymaker",
    hltMuonInputString_     = cms.string("unused"),
    hltElectronInputString_ = cms.string("IsoElectronProducer"),
    hltPfJetsInputTag_      = cms.InputTag("hltAK4PFJetsCorrected"),
    hltPfMetInputTag_       = cms.InputTag("hltPFMETProducer"),
    hltPfHTInputTag_        = cms.InputTag("hltPFHT"),
    hltCaloHTInputTag_      = cms.InputTag("hltHtMht"),
    hltGenJetsInputTag_     = cms.InputTag("ak4GenJetsNoNu"),
    recoMuonInputTag_       = cms.InputTag("unused"),
    recoElectronInputTag_   = cms.InputTag("unused"), 
    recoPfJetsInputTag_     = cms.InputTag("unused"),
    recoPfMetInputTag_      = cms.InputTag("unused"),
    recoGenJetsInputTag_    = cms.InputTag("unused"),
    m_Jets_                 = cms.InputTag("hltSelector4CentralJetsL1FastJet"),
    m_JetTags_              = cms.InputTag("hltCombinedSecondaryVertexBJetTagsCalo")
)
