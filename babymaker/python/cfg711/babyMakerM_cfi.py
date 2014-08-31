import FWCore.ParameterSet.Config as cms

babyMakerM = cms.EDProducer("BabyMakerM",
                            MuonInputTag_           = cms.InputTag("hltMuons"),
                            pfJetsInputTag_         = cms.InputTag("hltAK4PFJetsCorrected"),
                            pfMetInputTag_          = cms.InputTag("hltPFMETProducer"),
                            pfHTInputTag_           = cms.InputTag("hltPFHT"),
                            genJetsInputTag_        = cms.InputTag("ak4GenJetsNoNu"),
)
