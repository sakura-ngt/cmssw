import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer

# Tracker-alignment sensitive quantities from scouting di-muon resonances
# (charge-split mass bias, muon impact parameters by charge, PV-independent
# pair quantities, PV -> di-muon vertex displacement, pseudo-proper decay length).
# Histograms go to <FolderName>/DiMuonAlignmentMonitor.
ScoutingDiMuonAlignmentMonitor = DQMEDAnalyzer("ScoutingDiMuonAlignmentMonitor",
                                               FolderName = cms.string('HLT/ScoutingOffline/DiMuon'),
                                               decayMotherName = cms.string('Z'),
                                               muons = cms.InputTag('hltScoutingMuonPackerVtx'),
                                               primaryVertices = cms.InputTag('hltScoutingPrimaryVertexPacker', 'primaryVtx'),
                                               secondaryVertices = cms.InputTag('hltScoutingMuonPackerVtx', 'displacedVtx'),
                                               applyMuonID = cms.bool(True),
                                               minMuonPt = cms.double(3),
                                               maxMuonEta = cms.double(2.4),
                                               useClosestVertex = cms.bool(True),
                                               minVtxProb = cms.double(0.005),
                                               promptWindow = cms.double(200))

# J/psi: same observables, in its own folder; the pseudo-proper decay length is mostly meaningful here
ScoutingDiMuonAlignmentMonitorJpsi = ScoutingDiMuonAlignmentMonitor.clone(FolderName = 'HLT/ScoutingOffline/DiMuonJpsi',
                                                                           decayMotherName = 'J/#psi')
