import FWCore.ParameterSet.Config as cms

babyMakerE = cms.EDProducer("BabyMakerE",
                            ElectronInputTag_       = cms.InputTag("hltEgammaGsfElectrons"), 
                            pfJetsInputTag_         = cms.InputTag("hltAK4PFJetsCorrected"),
                            pfMetInputTag_          = cms.InputTag("hltPFMETProducer"),
                            pfHTInputTag_           = cms.InputTag("hltPFHT"),
                            genJetsInputTag_        = cms.InputTag("ak4GenJetsNoNu"),
)


