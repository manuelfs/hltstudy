import FWCore.ParameterSet.Config as cms

IsoMuonProducer = cms.EDProducer(
    "IsoMuonProducer",
    isoTag = cms.InputTag("hltL3MuonCombRelIsolationVVVL","combinedRelativeIsoDeposits"),
    candTag = cms.InputTag("hltL3MuonCandidates")
)
