// -*- C++ -*-
//
// Package:    HLTriggerOffline/Scouting
// Class:      ScoutingDiMuonVertexMonitor
//
// Description:
//   DQM monitor for di-muon vertices in Run 3 scouting data.
//   Inspired by DQMOffline/Alignment/src/DiMuonVertexMonitor.cc but adapted
//   for scouting objects:
//     - Uses Run3ScoutingMuon (no reco::Track, no TransientTrackBuilder)
//     - Uses Run3ScoutingVertex for both the pre-fitted di-muon secondary
//       vertex and the primary vertex collection
//     - Avoids KalmanVertexFitter (vertex already computed online)
//
//   For each opposite-sign muon pair inside the mass window of the decay
//   mother, the di-muon vertex fitted at HLT for that pair is taken from the
//   displaced-vertex collection (via Run3ScoutingMuon::vtxIndx()), and the
//   primary vertex closest to it is used as reference (or the first one if
//   useClosestVertex is false). Monitored quantities:
//     * Invariant mass of all opposite-sign pairs
//     * SV probability, chi2, chi2/ndof
//     * PV-SV distance (xy and 3D), its uncertainty (full covariance) and significance
//     * cos(phi) (xy and 3D) between the di-muon momentum and the PV->SV
//       displacement, for pairs with PV-SV xy distance below maxSVdistXY
//     * Per-muon track d_xy and d_z w.r.t. the selected primary vertex
//   Independently of the pair selection, the positions of all valid di-muon
//   vertices of the event are histogrammed in the lab frame (x vs y, and signed
//   r = sign(y) sqrt(x^2 + y^2) vs z), over the whole tracker and zoomed on the
//   beam pipe / pixel region.

// system includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

// CMSSW framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

// DQM
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

// Scouting data formats
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"

// Math
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TMath.h"

// Local utilities (ID selectors shared with ScoutingDileptonMonitor)
#include "ScoutingDQMUtils.h"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
namespace {
  inline math::PtEtaPhiMLorentzVector muonP4(const Run3ScoutingMuon& mu) {
    return math::PtEtaPhiMLorentzVector(mu.pt(), mu.eta(), mu.phi(), scoutingDQMUtils::MUON_MASS);
  }

  // Sum of the position covariance matrices of two vertices, (x, y, z) ordering.
  struct Cov3 {
    double xx, yy, zz, xy, xz, yz;
  };

  inline Cov3 sumCov(const Run3ScoutingVertex& a, const Run3ScoutingVertex& b) {
    return {a.xError() * a.xError() + b.xError() * b.xError(),
            a.yError() * a.yError() + b.yError() * b.yError(),
            a.zError() * a.zError() + b.zError() * b.zError(),
            static_cast<double>(a.xyCov()) + b.xyCov(),
            static_cast<double>(a.xzCov()) + b.xzCov(),
            static_cast<double>(a.yzCov()) + b.yzCov()};
  }

  // Transverse distance between two vertices and its uncertainty,
  // sigma^2 = d^T C d / |d|^2 with C the sum of the two covariance matrices.
  inline void vtxDistXY(const Run3ScoutingVertex& sv, const Run3ScoutingVertex& pv, double& dist, double& distErr) {
    const double dx = sv.x() - pv.x();
    const double dy = sv.y() - pv.y();
    const Cov3 c = sumCov(sv, pv);
    dist = std::sqrt(dx * dx + dy * dy);
    const double var = dx * dx * c.xx + dy * dy * c.yy + 2. * dx * dy * c.xy;
    distErr = (dist > 0. && var > 0.) ? std::sqrt(var) / dist : 0.;
  }

  // Same in 3D.
  inline void vtxDist3D(const Run3ScoutingVertex& sv, const Run3ScoutingVertex& pv, double& dist, double& distErr) {
    const double dx = sv.x() - pv.x();
    const double dy = sv.y() - pv.y();
    const double dz = sv.z() - pv.z();
    const Cov3 c = sumCov(sv, pv);
    dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    const double var =
        dx * dx * c.xx + dy * dy * c.yy + dz * dz * c.zz + 2. * (dx * dy * c.xy + dx * dz * c.xz + dy * dz * c.yz);
    distErr = (dist > 0. && var > 0.) ? std::sqrt(var) / dist : 0.;
  }

  // Track impact parameters w.r.t. a point, as reco::TrackBase::dxy(point) / dz(point).
  // Run3ScoutingMuon::trk_dxy() / trk_dz() are computed w.r.t. (0,0,0), so they are
  // recomputed here from the stored reference point and momentum direction.
  inline double trkDxy(const Run3ScoutingMuon& mu, const Run3ScoutingVertex& pv) {
    return -(mu.trk_vx() - pv.x()) * std::sin(mu.trk_phi()) + (mu.trk_vy() - pv.y()) * std::cos(mu.trk_phi());
  }

  inline double trkDz(const Run3ScoutingMuon& mu, const Run3ScoutingVertex& pv) {
    const double cosPhi = std::cos(mu.trk_phi());
    const double sinPhi = std::sin(mu.trk_phi());
    const double pzOverPt = std::sinh(mu.trk_eta());
    return (mu.trk_vz() - pv.z()) - ((mu.trk_vx() - pv.x()) * cosPhi + (mu.trk_vy() - pv.y()) * sinPhi) * pzOverPt;
  }

  // Unit-safe cm -> um conversion factor
  constexpr double cmToUm = 1.e4;
}  // namespace

// ---------------------------------------------------------------------------
// Monitor class
// ---------------------------------------------------------------------------
class ScoutingDiMuonVertexMonitor : public DQMEDAnalyzer {
public:
  explicit ScoutingDiMuonVertexMonitor(const edm::ParameterSet&);
  ~ScoutingDiMuonVertexMonitor() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(edm::Event const&, edm::EventSetup const&) override;

private:
  const Run3ScoutingVertex* findClosestVertex(const Run3ScoutingVertex& sv,
                                              const Run3ScoutingVertexCollection& vertices) const;

  // ---- configuration -------------------------------------------------------
  const std::string folderName_;
  const std::string decayMotherName_;

  // Mass window determined by decayMotherName_
  std::pair<double, double> massLimits_;

  // Muon selection
  const bool applyMuonID_;
  const double minPt_;
  const double maxEta_;

  // Vertex selection / quality cuts
  const bool useClosestVertex_;  // use the PV closest to the di-muon vertex instead of the first one
  const double minVtxProb_;
  const double maxSVdistXY_;  // um, cut on PV-SV xy distance for the pointing-angle histograms

  // Input tags / tokens
  const edm::EDGetTokenT<Run3ScoutingMuonCollection> muonToken_;
  const edm::EDGetTokenT<Run3ScoutingVertexCollection> pvToken_;
  const edm::EDGetTokenT<Run3ScoutingVertexCollection> svToken_;  // di-muon displaced vtx

  // ---- histograms ----------------------------------------------------------

  // Vertex quality
  MonitorElement* hSVProb_{nullptr};
  MonitorElement* hSVChi2_{nullptr};
  MonitorElement* hSVNormChi2_{nullptr};

  // PV-SV distance
  MonitorElement* hSVDistXY_{nullptr};
  MonitorElement* hSVDistXYErr_{nullptr};
  MonitorElement* hSVDistXYSig_{nullptr};
  MonitorElement* hSVDist3D_{nullptr};
  MonitorElement* hSVDist3DErr_{nullptr};
  MonitorElement* hSVDist3DSig_{nullptr};

  // Pointing angle
  MonitorElement* hCosPhi_{nullptr};    // cos phi in xy
  MonitorElement* hCosPhi3D_{nullptr};  // cos phi in 3D
  MonitorElement* hCosPhiInv_{nullptr};
  MonitorElement* hCosPhiInv3D_{nullptr};
  MonitorElement* hCosPhiUnbalance_{nullptr};
  MonitorElement* hCosPhi3DUnbalance_{nullptr};

  // Invariant mass
  MonitorElement* hInvMass_{nullptr};

  // Per-muon track impact parameters w.r.t. the primary vertex
  MonitorElement* hTrkDxy_{nullptr};
  MonitorElement* hTrkDz_{nullptr};
  MonitorElement* hTrkDxyErr_{nullptr};
  MonitorElement* hTrkDzErr_{nullptr};

  // Number of muons and vertices per event (diagnostic)
  MonitorElement* hNMuons_{nullptr};
  MonitorElement* hNSV_{nullptr};

  // Di-muon vertex positions (all valid vertices of the event)
  MonitorElement* hSVXvsY_{nullptr};
  MonitorElement* hSVRvsZ_{nullptr};
  MonitorElement* hSVXvsYZoom_{nullptr};
  MonitorElement* hSVRvsZZoom_{nullptr};
};

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
ScoutingDiMuonVertexMonitor::ScoutingDiMuonVertexMonitor(const edm::ParameterSet& iConfig)
    : folderName_(iConfig.getParameter<std::string>("FolderName")),
      decayMotherName_(iConfig.getParameter<std::string>("decayMotherName")),
      applyMuonID_(iConfig.getParameter<bool>("applyMuonID")),
      minPt_(iConfig.getParameter<double>("minMuonPt")),
      maxEta_(iConfig.getParameter<double>("maxMuonEta")),
      useClosestVertex_(iConfig.getParameter<bool>("useClosestVertex")),
      minVtxProb_(iConfig.getParameter<double>("minVtxProb")),
      maxSVdistXY_(iConfig.getParameter<double>("maxSVdistXY")),
      muonToken_(consumes<Run3ScoutingMuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      pvToken_(consumes<Run3ScoutingVertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      svToken_(consumes<Run3ScoutingVertexCollection>(iConfig.getParameter<edm::InputTag>("secondaryVertices"))) {
  // Determine mass window from decay mother name
  if (decayMotherName_.find('Z') != std::string::npos) {
    massLimits_ = {50., 120.};
  } else if (decayMotherName_.find("J/#psi") != std::string::npos ||
             decayMotherName_.find("J/psi") != std::string::npos) {
    massLimits_ = {2.7, 3.4};
  } else if (decayMotherName_.find("#Upsilon") != std::string::npos ||
             decayMotherName_.find("Upsilon") != std::string::npos) {
    massLimits_ = {8.9, 9.9};
  } else {
    edm::LogWarning("ScoutingDiMuonVertexMonitor")
        << "Unrecognised decay mother '" << decayMotherName_ << "'. Defaulting to Z window [50, 120] GeV.";
    massLimits_ = {50., 120.};
  }
}

// ---------------------------------------------------------------------------
// bookHistograms
// ---------------------------------------------------------------------------
void ScoutingDiMuonVertexMonitor::bookHistograms(DQMStore::IBooker& iBooker, edm::Run const&, edm::EventSetup const&) {
  iBooker.setCurrentFolder(folderName_ + "/DiMuonVertexMonitor");

  const std::string motName = decayMotherName_;
  const std::string ps = "N(#mu#mu pairs)";
  const std::string histTit = motName + " #rightarrow #mu^{+}#mu^{-}";

  // ---- vertex quality ----
  hSVProb_ = iBooker.book1D("VtxProb", ";" + motName + " vertex probability;" + ps, 100, 0., 1.);

  hSVChi2_ = iBooker.book1D("VtxChi2", ";#chi^{2} of the " + motName + " vertex;" + ps, 200, 0., 200.);

  hSVNormChi2_ = iBooker.book1D("VtxNormChi2", ";#chi^{2}/ndof of the " + motName + " vertex;" + ps, 100, 0., 20.);

  // ---- PV-SV distance xy ----
  hSVDistXY_ = iBooker.book1D("VtxDistXY", histTit + ";PV-" + motName + "V xy distance [#mum];" + ps, 100, 0., 300.);

  hSVDistXYErr_ =
      iBooker.book1D("VtxDistXYErr", histTit + ";PV-" + motName + "V xy distance error [#mum];" + ps, 100, 0., 1000.);

  hSVDistXYSig_ =
      iBooker.book1D("VtxDistXYSig", histTit + ";PV-" + motName + "V xy distance significance;" + ps, 100, 0., 20.);

  // ---- PV-SV distance 3D ----
  hSVDist3D_ = iBooker.book1D("VtxDist3D", histTit + ";PV-" + motName + "V 3D distance [#mum];" + ps, 100, 0., 1000.);

  hSVDist3DErr_ =
      iBooker.book1D("VtxDist3DErr", histTit + ";PV-" + motName + "V 3D distance error [#mum];" + ps, 100, 0., 1000.);

  hSVDist3DSig_ =
      iBooker.book1D("VtxDist3DSig", histTit + ";PV-" + motName + "V 3D distance significance;" + ps, 100, 0., 20.);

  // ---- pointing angle ----
  hCosPhi_ = iBooker.book1D("CosPhi", histTit + ";cos(#phi_{xy});" + ps, 50, -1., 1.);

  hCosPhi3D_ = iBooker.book1D("CosPhi3D", histTit + ";cos(#phi_{3D});" + ps, 50, -1., 1.);

  hCosPhiInv_ = iBooker.book1D("CosPhiInv", histTit + ";inverted cos(#phi_{xy});" + ps, 50, -1., 1.);

  hCosPhiInv3D_ = iBooker.book1D("CosPhiInv3D", histTit + ";inverted cos(#phi_{3D});" + ps, 50, -1., 1.);

  hCosPhiUnbalance_ = iBooker.book1D("CosPhiUnbalance", histTit + ";cos(#phi_{xy}) unbalance;#Delta" + ps, 50, -1., 1.);

  hCosPhi3DUnbalance_ =
      iBooker.book1D("CosPhi3DUnbalance", histTit + ";cos(#phi_{3D}) unbalance;#Delta" + ps, 50, -1., 1.);

  // ---- invariant mass ----
  hInvMass_ =
      iBooker.book1D("InvMass", histTit + ";M(#mu^{+}#mu^{-}) [GeV];" + ps, 70, massLimits_.first, massLimits_.second);

  // ---- track impact parameters ----
  hTrkDxy_ = iBooker.book1D("TrkDxy", histTit + ";muon track d_{xy}(PV) [#mum];muon tracks", 150, -300., 300.);

  hTrkDz_ = iBooker.book1D("TrkDz", histTit + ";muon track d_{z}(PV) [#mum];muon tracks", 150, -300., 300.);

  hTrkDxyErr_ = iBooker.book1D("TrkDxyErr", histTit + ";muon track err_{dxy} [#mum];muon tracks", 250, 0., 500.);

  hTrkDzErr_ = iBooker.book1D("TrkDzErr", histTit + ";muon track err_{dz} [#mum];muon tracks", 250, 0., 500.);

  // ---- diagnostics ----
  hNMuons_ = iBooker.book1D("NMuons", ";Number of selected muons;Events", 20, 0., 20.);
  hNSV_ = iBooker.book1D("NSV", ";Number of di-muon secondary vertices;Events", 20, 0., 20.);

  // ---- di-muon vertex positions (lab frame, cm) ----
  const std::string vtxTit = "all di-muon vertices";
  hSVXvsY_ = iBooker.book2D(
      "VtxXvsY", vtxTit + ";di-muon vertex x [cm];di-muon vertex y [cm]", 220, -110., 110., 220, -110., 110.);
  hSVRvsZ_ = iBooker.book2D("VtxRvsZ",
                            vtxTit + ";di-muon vertex z [cm];di-muon vertex signed r = sign(y)#sqrt{x^{2}+y^{2}} [cm]",
                            260,
                            -130.,
                            130.,
                            220,
                            -110.,
                            110.);
  hSVXvsYZoom_ = iBooker.book2D(
      "VtxXvsY_zoom", vtxTit + ";di-muon vertex x [cm];di-muon vertex y [cm]", 200, -20., 20., 200, -20., 20.);
  hSVRvsZZoom_ =
      iBooker.book2D("VtxRvsZ_zoom",
                     vtxTit + ";di-muon vertex z [cm];di-muon vertex signed r = sign(y)#sqrt{x^{2}+y^{2}} [cm]",
                     200,
                     -40.,
                     40.,
                     200,
                     -20.,
                     20.);
}

// ---------------------------------------------------------------------------
// Primary vertex closest (in 3D) to the di-muon vertex
// ---------------------------------------------------------------------------
const Run3ScoutingVertex* ScoutingDiMuonVertexMonitor::findClosestVertex(
    const Run3ScoutingVertex& sv, const Run3ScoutingVertexCollection& vertices) const {
  const Run3ScoutingVertex* closest = nullptr;
  double minDist2 = std::numeric_limits<double>::max();
  for (const auto& pv : vertices) {
    if (!pv.isValidVtx())
      continue;
    const double dx = sv.x() - pv.x();
    const double dy = sv.y() - pv.y();
    const double dz = sv.z() - pv.z();
    const double d2 = dx * dx + dy * dy + dz * dz;
    if (d2 < minDist2) {
      minDist2 = d2;
      closest = &pv;
    }
  }
  return closest;
}

// ---------------------------------------------------------------------------
// analyze
// ---------------------------------------------------------------------------
void ScoutingDiMuonVertexMonitor::analyze(edm::Event const& iEvent, edm::EventSetup const&) {
  // ---- Muon collection ----
  edm::Handle<Run3ScoutingMuonCollection> muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);
  if (!muonHandle.isValid()) {
    edm::LogWarning("ScoutingDiMuonVertexMonitor") << "Invalid muon collection handle.";
    return;
  }

  // ---- Primary vertex collection ----
  edm::Handle<Run3ScoutingVertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);
  if (!pvHandle.isValid() || pvHandle->empty()) {
    edm::LogWarning("ScoutingDiMuonVertexMonitor") << "No primary vertices found; skipping event.";
    return;
  }

  // ---- Secondary vertex collection (di-muon) ----
  edm::Handle<Run3ScoutingVertexCollection> svHandle;
  iEvent.getByToken(svToken_, svHandle);
  if (!svHandle.isValid()) {
    edm::LogWarning("ScoutingDiMuonVertexMonitor") << "Invalid secondary vertex collection handle.";
    return;
  }
  const auto& svs = *svHandle;

  hNSV_->Fill(static_cast<double>(svs.size()));

  // ---- Positions of all valid di-muon vertices (independent of the pair selection) ----
  for (const auto& sv : svs) {
    if (!sv.isValidVtx())
      continue;
    const double r = std::copysign(std::hypot(sv.x(), sv.y()), sv.y());
    hSVXvsY_->Fill(sv.x(), sv.y());
    hSVRvsZ_->Fill(sv.z(), r);
    hSVXvsYZoom_->Fill(sv.x(), sv.y());
    hSVRvsZZoom_->Fill(sv.z(), r);
  }

  // ---- Muon selection ----
  std::vector<size_t> selectedMuons;
  selectedMuons.reserve(muonHandle->size());
  for (size_t i = 0; i < muonHandle->size(); ++i) {
    const auto& mu = (*muonHandle)[i];
    if (mu.pt() < minPt_)
      continue;
    if (std::abs(mu.eta()) > maxEta_)
      continue;
    if (applyMuonID_ && !scoutingDQMUtils::scoutingMuonID(mu))
      continue;
    selectedMuons.push_back(i);
  }

  hNMuons_->Fill(static_cast<double>(selectedMuons.size()));

  if (selectedMuons.size() < 2)
    return;

  // muons whose impact parameters have already been filled in this event
  std::vector<bool> ipFilled(muonHandle->size(), false);

  // ---- Loop over opposite-sign muon pairs ----
  for (size_t ii = 0; ii < selectedMuons.size(); ++ii) {
    for (size_t jj = ii + 1; jj < selectedMuons.size(); ++jj) {
      const auto& mu_i = (*muonHandle)[selectedMuons[ii]];
      const auto& mu_j = (*muonHandle)[selectedMuons[jj]];

      // Opposite sign
      if (mu_i.charge() * mu_j.charge() >= 0)
        continue;

      // Invariant mass (all opposite-sign pairs), then the mass window of the decay mother
      const auto p4pair = muonP4(mu_i) + muonP4(mu_j);
      const double mass = p4pair.mass();
      hInvMass_->Fill(mass);
      if (mass < massLimits_.first || mass > massLimits_.second)
        continue;

      // ---- Di-muon vertex fitted at HLT for this pair ----
      // Run3ScoutingMuon::vtxIndx() holds the indices, in the displaced-vertex
      // collection, of the vertices the muon belongs to: the pair's vertex is a
      // common index of the two muons (the most probable one if there are several).
      const Run3ScoutingVertex* bestSV = nullptr;
      double bestProb = -1.;
      for (const int idx : mu_i.vtxIndx()) {
        if (idx < 0 || static_cast<size_t>(idx) >= svs.size())
          continue;
        const auto& idxJ = mu_j.vtxIndx();
        if (std::find(idxJ.begin(), idxJ.end(), idx) == idxJ.end())
          continue;
        const auto& sv = svs[idx];
        if (!sv.isValidVtx())
          continue;
        const double prob = (sv.ndof() > 0) ? TMath::Prob(sv.chi2(), sv.ndof()) : 0.;
        if (prob > bestProb) {
          bestProb = prob;
          bestSV = &sv;
        }
      }

      if (bestSV == nullptr)
        continue;  // no vertex fitted for this pair

      // ---- Vertex probability and chi2 ----
      const double chi2 = bestSV->chi2();
      const int ndof = bestSV->ndof();

      hSVProb_->Fill(bestProb);
      hSVChi2_->Fill(chi2);
      if (ndof > 0)
        hSVNormChi2_->Fill(chi2 / ndof);

      if (bestProb < minVtxProb_)
        continue;

      // ---- Reference primary vertex ----
      const Run3ScoutingVertex* closestPV = useClosestVertex_ ? findClosestVertex(*bestSV, *pvHandle) : nullptr;
      const Run3ScoutingVertex& pv = closestPV ? *closestPV : pvHandle->front();

      // ---- Per-muon impact parameters w.r.t. the primary vertex (once per muon) ----
      for (const size_t k : {selectedMuons[ii], selectedMuons[jj]}) {
        if (ipFilled[k])
          continue;
        ipFilled[k] = true;
        const auto& mu = (*muonHandle)[k];
        hTrkDxy_->Fill(trkDxy(mu, pv) * cmToUm);
        hTrkDz_->Fill(trkDz(mu, pv) * cmToUm);
        hTrkDxyErr_->Fill(mu.trk_dxyError() * cmToUm);
        hTrkDzErr_->Fill(mu.trk_dzError() * cmToUm);
      }

      // ---- PV-SV distances (full covariance of both vertices) ----
      double distXY, distXYErr;
      vtxDistXY(*bestSV, pv, distXY, distXYErr);

      double dist3D, dist3DErr;
      vtxDist3D(*bestSV, pv, dist3D, dist3DErr);

      hSVDistXY_->Fill(distXY * cmToUm);
      hSVDistXYErr_->Fill(distXYErr * cmToUm);
      if (distXYErr > 0.)
        hSVDistXYSig_->Fill(distXY / distXYErr);

      hSVDist3D_->Fill(dist3D * cmToUm);
      hSVDist3DErr_->Fill(dist3DErr * cmToUm);
      if (dist3DErr > 0.)
        hSVDist3DSig_->Fill(dist3D / dist3DErr);

      // ---- Pointing angles, for vertices close to the PV (as in DiMuonVertexMonitor) ----
      if (distXY * cmToUm > maxSVdistXY_)
        continue;

      // Displacement vector PV -> SV
      const double dVtxX = bestSV->x() - pv.x();
      const double dVtxY = bestSV->y() - pv.y();
      const double dVtxZ = bestSV->z() - pv.z();

      const double pxDimu = p4pair.px();
      const double pyDimu = p4pair.py();
      const double pzDimu = p4pair.pz();

      const double magPtDimu = std::sqrt(pxDimu * pxDimu + pyDimu * pyDimu);
      const double magDVtxXY = std::sqrt(dVtxX * dVtxX + dVtxY * dVtxY);
      const double magPDimu = std::sqrt(pxDimu * pxDimu + pyDimu * pyDimu + pzDimu * pzDimu);
      const double magDVtx3D = std::sqrt(dVtxX * dVtxX + dVtxY * dVtxY + dVtxZ * dVtxZ);

      // cos phi in the xy plane
      if (magPtDimu > 0. && magDVtxXY > 0.) {
        const double cosphi = (pxDimu * dVtxX + pyDimu * dVtxY) / (magPtDimu * magDVtxXY);
        hCosPhi_->Fill(cosphi);
        hCosPhiInv_->Fill(-cosphi);
        hCosPhiUnbalance_->Fill(cosphi, 1.);
        hCosPhiUnbalance_->Fill(-cosphi, -1.);
      }

      // cos phi in 3D
      if (magPDimu > 0. && magDVtx3D > 0.) {
        const double cosphi3D = (pxDimu * dVtxX + pyDimu * dVtxY + pzDimu * dVtxZ) / (magPDimu * magDVtx3D);
        hCosPhi3D_->Fill(cosphi3D);
        hCosPhiInv3D_->Fill(-cosphi3D);
        hCosPhi3DUnbalance_->Fill(cosphi3D, 1.);
        hCosPhi3DUnbalance_->Fill(-cosphi3D, -1.);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// fillDescriptions
// ---------------------------------------------------------------------------
void ScoutingDiMuonVertexMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("FolderName", "HLT/ScoutingOffline/DiMuon");
  desc.add<std::string>("decayMotherName", "Z");

  desc.add<edm::InputTag>("muons", edm::InputTag("hltScoutingMuonPackerVtx"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"));
  desc.add<edm::InputTag>("secondaryVertices", edm::InputTag("hltScoutingMuonPackerVtx", "displacedVtx"));

  desc.add<bool>("applyMuonID", true);
  desc.add<double>("minMuonPt", 3.0);
  desc.add<double>("maxMuonEta", 2.4);
  desc.add<bool>("useClosestVertex", true)
      ->setComment("use the primary vertex closest to the di-muon vertex instead of the first one");
  desc.add<double>("minVtxProb", 0.005);
  desc.add<double>("maxSVdistXY", 50.0)->setComment("um, max PV-SV xy distance for the pointing-angle histograms");

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(ScoutingDiMuonVertexMonitor);
