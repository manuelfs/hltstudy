import FWCore.ParameterSet.Config as cms

IsoElectronProducer = cms.EDProducer(
    "IsoElectronProducer",
    ecalTag = cms.InputTag("hltEle27WP80EcalIsoFilter"),
    hcalTag = cms.InputTag("hltEle27WP80HcalIsoFilter"),
    trackTag = cms.InputTag("hltEle27WP80GsfTrackIsoFilter"),
    ecalMapTag = cms.InputTag("hltEgammaEcalPFClusterIso"),
    hcalMapTag = cms.InputTag("hltEgammaHcalPFClusterIso"),
    trackMapTag = cms.InputTag("hltEgammaEleGsfTrackIso")
)
