import FWCore.ParameterSet.Config as cms

IsoElectronProducer = cms.EDProducer(
    "IsoElectronProducer",
    clustershapeTag = cms.InputTag("hltEle15VVVLClusterShapeFilter"),
    clustershapeMapTag = cms.InputTag("hltEgammaClusterShape"),
    heTag = cms.InputTag("hltEle15VVVLHEFilter"),
    heMapTag = cms.InputTag("hltEgammaHoverE"),
    eminuspTag = cms.InputTag("hltEle15VVVLGsfOneOEMinusOneOPFilter"),
    eminuspMapTag = cms.InputTag('hltEgammaGsfTrackVars','OneOESuperMinusOneOP'),
    detaTag = cms.InputTag("hltEle15VVVLGsfDetaFilter"),
    detaMapTag = cms.InputTag('hltEgammaGsfTrackVars','Deta'),
    dphiTag = cms.InputTag("hltEle15VVVLGsfDphiFilter"),
    dphiMapTag = cms.InputTag('hltEgammaGsfTrackVars','Dphi'),

    ecalTag = cms.InputTag("hltEle15VVVLEcalIsoFilter"),
    hcalTag = cms.InputTag("hltEle15VVVLHcalIsoFilter"),
    trackTag = cms.InputTag("hltEle15VVVLGsfTrackIsoFilter"),
    ecalMapTag = cms.InputTag("hltEgammaEcalPFClusterIso"),
    hcalMapTag = cms.InputTag("hltEgammaHcalPFClusterIso"),
    trackMapTag = cms.InputTag("hltEgammaEleGsfTrackIso")
)
