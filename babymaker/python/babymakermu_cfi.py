import FWCore.ParameterSet.Config as cms

babymaker = cms.EDProducer(
    "babymaker",
    hltMuonInputTag_       = cms.InputTag("hltMuons"),
    hltElectronInputTag_   = cms.InputTag("unused"), 
    hltPfJetsInputTag_     = cms.InputTag("hltAK4PFJetsCorrected"),
    #hltPfJetsInputTag_     = cms.InputTag("hltAK4PFJetL1FastL2L3CorrectedNoPU"),
    hltPfMetInputTag_      = cms.InputTag("hltPFMETProducer"),
    hltPfHTInputTag_       = cms.InputTag("hltPFHT"),
    #hltPfHTInputTag_       = cms.InputTag("hltPFHTNoPU"),
    hltGenJetsInputTag_    = cms.InputTag("ak4GenJetsNoNu"),
    recoMuonInputTag_      = cms.InputTag("muons"),
    recoElectronInputTag_  = cms.InputTag("gedGsfElectrons"), 
    recoPfJetsInputTag_    = cms.InputTag("ak4PFJets"),
    recoPfMetInputTag_     = cms.InputTag("pfMet"),
    recoGenJetsInputTag_   = cms.InputTag("ak5GenJets"),
    m_Jets_                = cms.InputTag("hltSelector4CentralJetsL1FastJet"), # btagging 
    m_JetTags_             = cms.InputTag("hltL3CombinedSecondaryVertexBJetTags"), # btagging 
)
