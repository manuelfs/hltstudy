import FWCore.ParameterSet.Config as cms

babymaker = cms.EDProducer(
    "babymaker",
    hltMuonInputString_     = cms.string("unused"),
    hltElectronInputString_ = cms.string("IsoElectronProducer"),
    hltPfJetsInputTag_      = cms.InputTag("hltAK4PFJetsCorrected"),
    hltPfMetInputTag_       = cms.InputTag("hltPFMETProducer"),
    hltPfHTInputTag_        = cms.InputTag("hltPFHT"),
    hltCaloMetInputTag_     = cms.InputTag("hltHtMht"),
    hltCaloHTInputTag_      = cms.InputTag("hltHtMht"),
    hltGenJetsInputTag_     = cms.InputTag("ak4GenJetsNoNu"),
    recoMuonInputTag_       = cms.InputTag("muons"),
    recoElectronInputTag_   = cms.InputTag("gedGsfElectrons"), 
    recoPfJetsInputTag_     = cms.InputTag("ak4PFJets"),
    recoPfMetInputTag_      = cms.InputTag("pfMet"),
    recoGenJetsInputTag_    = cms.InputTag("ak5GenJets"),
    m_Jets_                 = cms.InputTag("hltL3CombinedSecondaryVertexBJetTags"),
    m_JetTags_              = cms.InputTag("hltL3CombinedSecondaryVertexBJetTags")
)
