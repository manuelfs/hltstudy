import FWCore.ParameterSet.Config as cms

IsoMuonProducer = cms.EDProducer(
    "IsoMuonProducer",
    isoTag = cms.InputTag("hltL3MuonCombRelIsolationsIterTrkRegIter02","combinedRelativeIsoDeposits"),
    candTag = cms.InputTag("hltL3MuonCandidates")
)
