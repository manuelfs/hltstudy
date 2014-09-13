import FWCore.ParameterSet.Config as cms

IsoElectronProducer = cms.EDProducer(
    "IsoElectronProducer",
    ecalTag = cms.InputTag("hltEle15IsoVVVLEcalIsoFilter"),
    hcalTag = cms.InputTag("hltEle15IsoVVVLHcalIsoFilter"),
    trackTag = cms.InputTag("hltEle15IsoVVVLGsfTrackIsoFilter"),
    ecalMapTag = cms.InputTag("hltEgammaEcalPFClusterIso"),
    hcalMapTag = cms.InputTag("hltEgammaHcalPFClusterIso"),
    trackMapTag = cms.InputTag("hltEgammaEleGsfTrackIso")
)
