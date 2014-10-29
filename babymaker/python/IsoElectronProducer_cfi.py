import FWCore.ParameterSet.Config as cms

IsoElectronProducer = cms.EDProducer(
    "IsoElectronProducer",
    clustershapeTag = cms.InputTag("hltEle15IsoVVVLClusterShapeFilter"),
    clustershapeMapTag = cms.InputTag("hltEgammaClusterShape"),
    heTag = cms.InputTag("hltEle15IsoVVVLHEFilter"),
    heMapTag = cms.InputTag("hltEgammaHoverE"),
    eminuspTag = cms.InputTag("hltEle15IsoVVVLGsfOneOEMinusOneOPFilter"),
    eminuspMapTag = cms.InputTag('hltEgammaGsfTrackVars','OneOESuperMinusOneOP'),
    detaTag = cms.InputTag("hltEle15IsoVVVLGsfDetaFilter"),
    detaMapTag = cms.InputTag('hltEgammaGsfTrackVars','Deta'),
    dphiTag = cms.InputTag("hltEle15IsoVVVLGsfDphiFilter"),
    dphiMapTag = cms.InputTag('hltEgammaGsfTrackVars','Dphi'),

    ecalTag = cms.InputTag("hltEle15IsoVVVLEcalIsoFilter"),
    hcalTag = cms.InputTag("hltEle15IsoVVVLHcalIsoFilter"),
    trackTag = cms.InputTag("hltEle15IsoVVVLGsfTrackIsoFilter"),
    ecalMapTag = cms.InputTag("hltEgammaEcalPFClusterIso"),
    hcalMapTag = cms.InputTag("hltEgammaHcalPFClusterIso"),
    trackMapTag = cms.InputTag("hltEgammaEleGsfTrackIso")
)
