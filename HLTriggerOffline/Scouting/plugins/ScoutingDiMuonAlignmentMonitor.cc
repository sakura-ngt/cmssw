// -*- C++ -*-
//
// Package:    HLTriggerOffline/Scouting
// Class:      ScoutingDiMuonAlignmentMonitor
//
// Description:
//   DQM monitor of tracker-alignment sensitive quantities built from di-muon
//   resonances (Z, J/psi, Upsilon) in Run 3 scouting data. It complements
//   ScoutingDiMuonVertexMonitor with observables aimed at the tracker weak
//   modes (deformations DeltaR, Deltaphi, Deltaz as a function of R, z, phi):
//
//   1. Charge-split mass bias (Deltaphi modes: layer rotation, twist, sagitta,
//      i.e. charge-antisymmetric curvature biases):
//        di-muon mass vs eta and phi of mu+ and of mu-, vs eta(mu+) - eta(mu-),
//        and vs (eta(mu+), eta(mu-)); as 2D histograms (mass distribution per bin,
//        for fits in harvesting / offline) and as profiles (quick look).
//   2. Muon impact parameters w.r.t. the primary vertex, split by charge:
//        <d_xy> vs phi and <d_z> vs eta (elliptical, skew, telescope, bowing).
//   3. Primary-vertex independent pair quantities: differences of the d_xy / d_z
//      of the two muons, both measured w.r.t. their common di-muon vertex,
//      vs eta / phi (relative misalignment between detector regions).
//   4. PV -> di-muon vertex displacement components vs eta / phi of the pair:
//        Delta z, L_xy (along the pair p_T) and L_perp (across it) (Deltaz modes,
//        beam-spot and longitudinal biases).
//   5. Pseudo-proper decay length l_xy = L_xy * m / p_T, its uncertainty (full
//      vertex covariance) and significance, and <l_xy> of the prompt peak vs
//      eta / phi of the pair (mainly for J/psi: the prompt peak is centred at 0
//      for a well aligned tracker).
//
//   Pairs: opposite-sign muons passing the muon selection, inside the mass window
//   of the decay mother. The mass-bias histograms use all such pairs; the vertex
//   based ones (2-5) require the di-muon vertex fitted at HLT for the pair (from
//   Run3ScoutingMuon::vtxIndx()) with probability above minVtxProb, and use the
//   primary vertex closest to it (or the first one if useClosestVertex is false).

// system includes
#include <algorithm>
#include <cmath>
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

// Local utilities (muon ID shared with the other scouting monitors)
#include "ScoutingDQMUtils.h"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
namespace {
  constexpr double cmToUm = 1.e4;

  inline math::PtEtaPhiMLorentzVector muonP4(const Run3ScoutingMuon& mu) {
    return math::PtEtaPhiMLorentzVector(mu.pt(), mu.eta(), mu.phi(), scoutingDQMUtils::MUON_MASS);
  }

  // Track impact parameters w.r.t. a point (x, y, z), as reco::TrackBase::dxy(point) / dz(point),
  // from the stored reference point and momentum direction of the muon track.
  inline double trkDxy(const Run3ScoutingMuon& mu, double x, double y) {
    return -(mu.trk_vx() - x) * std::sin(mu.trk_phi()) + (mu.trk_vy() - y) * std::cos(mu.trk_phi());
  }

  inline double trkDz(const Run3ScoutingMuon& mu, double x, double y, double z) {
    return (mu.trk_vz() - z) -
           ((mu.trk_vx() - x) * std::cos(mu.trk_phi()) + (mu.trk_vy() - y) * std::sin(mu.trk_phi())) *
               std::sinh(mu.trk_eta());
  }

  // Transverse position covariance (xx, yy, xy) of the sum of two vertices.
  struct Cov2 {
    double xx, yy, xy;
  };

  inline Cov2 sumCovXY(const Run3ScoutingVertex& a, const Run3ScoutingVertex& b) {
    return {a.xError() * a.xError() + b.xError() * b.xError(),
            a.yError() * a.yError() + b.yError() * b.yError(),
            static_cast<double>(a.xyCov()) + b.xyCov()};
  }
}  // namespace

// ---------------------------------------------------------------------------
// Monitor class
// ---------------------------------------------------------------------------
class ScoutingDiMuonAlignmentMonitor : public DQMEDAnalyzer {
public:
  explicit ScoutingDiMuonAlignmentMonitor(const edm::ParameterSet&);
  ~ScoutingDiMuonAlignmentMonitor() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(edm::Event const&, edm::EventSetup const&) override;

private:
  const Run3ScoutingVertex* findPairVertex(const Run3ScoutingMuon& mu1,
                                           const Run3ScoutingMuon& mu2,
                                           const Run3ScoutingVertexCollection& svs,
                                           double& prob) const;
  const Run3ScoutingVertex* findClosestVertex(const Run3ScoutingVertex& sv,
                                              const Run3ScoutingVertexCollection& pvs) const;

  // ---- configuration -------------------------------------------------------
  const std::string folderName_;
  const std::string decayMotherName_;
  std::pair<double, double> massLimits_;

  const bool applyMuonID_;
  const double minPt_;
  const double maxEta_;
  const bool useClosestVertex_;
  const double minVtxProb_;
  const double promptWindow_;  // um, |l_xy| window for the prompt-peak profiles

  const edm::EDGetTokenT<Run3ScoutingMuonCollection> muonToken_;
  const edm::EDGetTokenT<Run3ScoutingVertexCollection> pvToken_;
  const edm::EDGetTokenT<Run3ScoutingVertexCollection> svToken_;

  // ---- histograms ----------------------------------------------------------
  // 1. mass bias
  MonitorElement* hMassVsEtaPlus_{nullptr};
  MonitorElement* hMassVsEtaMinus_{nullptr};
  MonitorElement* hMassVsPhiPlus_{nullptr};
  MonitorElement* hMassVsPhiMinus_{nullptr};
  MonitorElement* hMassVsDeltaEta_{nullptr};
  MonitorElement* pMassVsEtaPlus_{nullptr};
  MonitorElement* pMassVsEtaMinus_{nullptr};
  MonitorElement* pMassVsPhiPlus_{nullptr};
  MonitorElement* pMassVsPhiMinus_{nullptr};
  MonitorElement* pMassVsDeltaEta_{nullptr};
  MonitorElement* pMassVsEtaPlusEtaMinus_{nullptr};

  // 2. impact parameters w.r.t. the PV, by charge
  MonitorElement* pDxyVsPhiPlus_{nullptr};
  MonitorElement* pDxyVsPhiMinus_{nullptr};
  MonitorElement* pDzVsEtaPlus_{nullptr};
  MonitorElement* pDzVsEtaMinus_{nullptr};

  // 3. PV-independent pair quantities (w.r.t. the di-muon vertex)
  MonitorElement* hDeltaDxy_{nullptr};
  MonitorElement* hDeltaDz_{nullptr};
  MonitorElement* pDeltaDxyVsPhiPlus_{nullptr};
  MonitorElement* pDeltaDzVsEtaPlus_{nullptr};
  MonitorElement* pDeltaDzVsEtaPlusEtaMinus_{nullptr};

  // 4. PV -> di-muon vertex displacement
  MonitorElement* pVtxDzVsEta_{nullptr};
  MonitorElement* pVtxDzVsPhi_{nullptr};
  MonitorElement* pLxyVsEta_{nullptr};
  MonitorElement* pLxyVsPhi_{nullptr};
  MonitorElement* pLperpVsPhi_{nullptr};

  // 5. pseudo-proper decay length
  MonitorElement* hDecayLength_{nullptr};
  MonitorElement* hDecayLengthErr_{nullptr};
  MonitorElement* hDecayLengthSig_{nullptr};
  MonitorElement* pDecayLengthVsEta_{nullptr};
  MonitorElement* pDecayLengthVsPhi_{nullptr};
};

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
ScoutingDiMuonAlignmentMonitor::ScoutingDiMuonAlignmentMonitor(const edm::ParameterSet& iConfig)
    : folderName_(iConfig.getParameter<std::string>("FolderName")),
      decayMotherName_(iConfig.getParameter<std::string>("decayMotherName")),
      applyMuonID_(iConfig.getParameter<bool>("applyMuonID")),
      minPt_(iConfig.getParameter<double>("minMuonPt")),
      maxEta_(iConfig.getParameter<double>("maxMuonEta")),
      useClosestVertex_(iConfig.getParameter<bool>("useClosestVertex")),
      minVtxProb_(iConfig.getParameter<double>("minVtxProb")),
      promptWindow_(iConfig.getParameter<double>("promptWindow")),
      muonToken_(consumes<Run3ScoutingMuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      pvToken_(consumes<Run3ScoutingVertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      svToken_(consumes<Run3ScoutingVertexCollection>(iConfig.getParameter<edm::InputTag>("secondaryVertices"))) {
  // same mass windows as ScoutingDiMuonVertexMonitor
  if (decayMotherName_.find('Z') != std::string::npos) {
    massLimits_ = {50., 120.};
  } else if (decayMotherName_.find("J/#psi") != std::string::npos ||
             decayMotherName_.find("J/psi") != std::string::npos) {
    massLimits_ = {2.7, 3.4};
  } else if (decayMotherName_.find("#Upsilon") != std::string::npos ||
             decayMotherName_.find("Upsilon") != std::string::npos) {
    massLimits_ = {8.9, 9.9};
  } else {
    edm::LogWarning("ScoutingDiMuonAlignmentMonitor")
        << "Unrecognised decay mother '" << decayMotherName_ << "'. Defaulting to Z window [50, 120] GeV.";
    massLimits_ = {50., 120.};
  }
}

// ---------------------------------------------------------------------------
// bookHistograms
// ---------------------------------------------------------------------------
void ScoutingDiMuonAlignmentMonitor::bookHistograms(DQMStore::IBooker& iBooker,
                                                    edm::Run const&,
                                                    edm::EventSetup const&) {
  iBooker.setCurrentFolder(folderName_ + "/DiMuonAlignmentMonitor");

  const std::string tit = decayMotherName_ + " #rightarrow #mu^{+}#mu^{-}";
  const std::string mLab = "M(#mu^{+}#mu^{-}) [GeV]";
  const double mMin = massLimits_.first;
  const double mMax = massLimits_.second;
  constexpr int nMass = 60;
  constexpr int nEta = 24;
  constexpr int nPhi = 24;
  constexpr double etaMax = 2.4;
  constexpr int nEta2 = 12;  // coarser binning for the (eta+, eta-) maps
  constexpr int nEtaPair = 24;
  constexpr double etaPairMax = 4.8;

  // ---- 1. mass bias ----
  hMassVsEtaPlus_ =
      iBooker.book2D("MassVsEtaPlus", tit + ";#mu^{+} #eta;" + mLab, nEta, -etaMax, etaMax, nMass, mMin, mMax);
  hMassVsEtaMinus_ =
      iBooker.book2D("MassVsEtaMinus", tit + ";#mu^{-} #eta;" + mLab, nEta, -etaMax, etaMax, nMass, mMin, mMax);
  hMassVsPhiPlus_ =
      iBooker.book2D("MassVsPhiPlus", tit + ";#mu^{+} #phi;" + mLab, nPhi, -M_PI, M_PI, nMass, mMin, mMax);
  hMassVsPhiMinus_ =
      iBooker.book2D("MassVsPhiMinus", tit + ";#mu^{-} #phi;" + mLab, nPhi, -M_PI, M_PI, nMass, mMin, mMax);
  hMassVsDeltaEta_ = iBooker.book2D("MassVsDeltaEta",
                                    tit + ";#eta(#mu^{+}) - #eta(#mu^{-});" + mLab,
                                    nEta,
                                    -2. * etaMax,
                                    2. * etaMax,
                                    nMass,
                                    mMin,
                                    mMax);

  pMassVsEtaPlus_ = iBooker.bookProfile(
      "p_MassVsEtaPlus", tit + ";#mu^{+} #eta;<" + mLab + ">", nEta, -etaMax, etaMax, mMin, mMax, "");
  pMassVsEtaMinus_ = iBooker.bookProfile(
      "p_MassVsEtaMinus", tit + ";#mu^{-} #eta;<" + mLab + ">", nEta, -etaMax, etaMax, mMin, mMax, "");
  pMassVsPhiPlus_ =
      iBooker.bookProfile("p_MassVsPhiPlus", tit + ";#mu^{+} #phi;<" + mLab + ">", nPhi, -M_PI, M_PI, mMin, mMax, "");
  pMassVsPhiMinus_ =
      iBooker.bookProfile("p_MassVsPhiMinus", tit + ";#mu^{-} #phi;<" + mLab + ">", nPhi, -M_PI, M_PI, mMin, mMax, "");
  pMassVsDeltaEta_ = iBooker.bookProfile("p_MassVsDeltaEta",
                                         tit + ";#eta(#mu^{+}) - #eta(#mu^{-});<" + mLab + ">",
                                         nEta,
                                         -2. * etaMax,
                                         2. * etaMax,
                                         mMin,
                                         mMax,
                                         "");
  pMassVsEtaPlusEtaMinus_ = iBooker.bookProfile2D("p_MassVsEtaPlusEtaMinus",
                                                  tit + ";#mu^{+} #eta;#mu^{-} #eta;<" + mLab + ">",
                                                  nEta2,
                                                  -etaMax,
                                                  etaMax,
                                                  nEta2,
                                                  -etaMax,
                                                  etaMax,
                                                  mMin,
                                                  mMax,
                                                  "");

  // ---- 2. impact parameters w.r.t. the PV, by charge ----
  constexpr double ipMax = 1000.;  // um, profile range
  pDxyVsPhiPlus_ = iBooker.bookProfile(
      "p_DxyVsPhiPlus", tit + ";#mu^{+} #phi;<d_{xy}(PV)> [#mum]", nPhi, -M_PI, M_PI, -ipMax, ipMax, "");
  pDxyVsPhiMinus_ = iBooker.bookProfile(
      "p_DxyVsPhiMinus", tit + ";#mu^{-} #phi;<d_{xy}(PV)> [#mum]", nPhi, -M_PI, M_PI, -ipMax, ipMax, "");
  pDzVsEtaPlus_ = iBooker.bookProfile(
      "p_DzVsEtaPlus", tit + ";#mu^{+} #eta;<d_{z}(PV)> [#mum]", nEta, -etaMax, etaMax, -ipMax, ipMax, "");
  pDzVsEtaMinus_ = iBooker.bookProfile(
      "p_DzVsEtaMinus", tit + ";#mu^{-} #eta;<d_{z}(PV)> [#mum]", nEta, -etaMax, etaMax, -ipMax, ipMax, "");

  // ---- 3. PV-independent pair quantities ----
  hDeltaDxy_ = iBooker.book1D("DeltaDxy",
                              tit + ";d_{xy}(#mu^{+}) - d_{xy}(#mu^{-}) w.r.t. #mu#mu vertex [#mum];N(#mu#mu pairs)",
                              100,
                              -200.,
                              200.);
  hDeltaDz_ = iBooker.book1D(
      "DeltaDz", tit + ";d_{z}(#mu^{+}) - d_{z}(#mu^{-}) w.r.t. #mu#mu vertex [#mum];N(#mu#mu pairs)", 100, -500., 500.);
  pDeltaDxyVsPhiPlus_ = iBooker.bookProfile("p_DeltaDxyVsPhiPlus",
                                            tit + ";#mu^{+} #phi;<d_{xy}(#mu^{+}) - d_{xy}(#mu^{-})> [#mum]",
                                            nPhi,
                                            -M_PI,
                                            M_PI,
                                            -ipMax,
                                            ipMax,
                                            "");
  pDeltaDzVsEtaPlus_ = iBooker.bookProfile("p_DeltaDzVsEtaPlus",
                                           tit + ";#mu^{+} #eta;<d_{z}(#mu^{+}) - d_{z}(#mu^{-})> [#mum]",
                                           nEta,
                                           -etaMax,
                                           etaMax,
                                           -ipMax,
                                           ipMax,
                                           "");
  pDeltaDzVsEtaPlusEtaMinus_ =
      iBooker.bookProfile2D("p_DeltaDzVsEtaPlusEtaMinus",
                            tit + ";#mu^{+} #eta;#mu^{-} #eta;<d_{z}(#mu^{+}) - d_{z}(#mu^{-})> [#mum]",
                            nEta2,
                            -etaMax,
                            etaMax,
                            nEta2,
                            -etaMax,
                            etaMax,
                            -ipMax,
                            ipMax,
                            "");

  // ---- 4. PV -> di-muon vertex displacement ----
  constexpr double dMax = 2000.;  // um, profile range
  pVtxDzVsEta_ = iBooker.bookProfile("p_VtxDzVsEta",
                                     tit + ";#mu#mu #eta;<z(#mu#mu vtx) - z(PV)> [#mum]",
                                     nEtaPair,
                                     -etaPairMax,
                                     etaPairMax,
                                     -dMax,
                                     dMax,
                                     "");
  pVtxDzVsPhi_ = iBooker.bookProfile(
      "p_VtxDzVsPhi", tit + ";#mu#mu #phi;<z(#mu#mu vtx) - z(PV)> [#mum]", nPhi, -M_PI, M_PI, -dMax, dMax, "");
  pLxyVsEta_ = iBooker.bookProfile(
      "p_LxyVsEta", tit + ";#mu#mu #eta;<L_{xy}> [#mum]", nEtaPair, -etaPairMax, etaPairMax, -dMax, dMax, "");
  pLxyVsPhi_ =
      iBooker.bookProfile("p_LxyVsPhi", tit + ";#mu#mu #phi;<L_{xy}> [#mum]", nPhi, -M_PI, M_PI, -dMax, dMax, "");
  pLperpVsPhi_ =
      iBooker.bookProfile("p_LperpVsPhi", tit + ";#mu#mu #phi;<L_{#perp}> [#mum]", nPhi, -M_PI, M_PI, -dMax, dMax, "");

  // ---- 5. pseudo-proper decay length ----
  hDecayLength_ =
      iBooker.book1D("DecayLength", tit + ";l_{xy} = L_{xy} m/p_{T} [#mum];N(#mu#mu pairs)", 200, -500., 1500.);
  hDecayLengthErr_ = iBooker.book1D("DecayLengthErr", tit + ";#sigma(l_{xy}) [#mum];N(#mu#mu pairs)", 100, 0., 200.);
  hDecayLengthSig_ = iBooker.book1D("DecayLengthSig", tit + ";l_{xy}/#sigma(l_{xy});N(#mu#mu pairs)", 120, -10., 50.);
  pDecayLengthVsEta_ = iBooker.bookProfile("p_DecayLengthVsEta",
                                           tit + ";#mu#mu #eta;<l_{xy}> prompt window [#mum]",
                                           nEtaPair,
                                           -etaPairMax,
                                           etaPairMax,
                                           -promptWindow_,
                                           promptWindow_,
                                           "");
  pDecayLengthVsPhi_ = iBooker.bookProfile("p_DecayLengthVsPhi",
                                           tit + ";#mu#mu #phi;<l_{xy}> prompt window [#mum]",
                                           nPhi,
                                           -M_PI,
                                           M_PI,
                                           -promptWindow_,
                                           promptWindow_,
                                           "");
}

// ---------------------------------------------------------------------------
// Di-muon vertex fitted at HLT for a pair: common index of the two muons in
// the displaced-vertex collection (Run3ScoutingMuon::vtxIndx()), the most
// probable one if there are several.
// ---------------------------------------------------------------------------
const Run3ScoutingVertex* ScoutingDiMuonAlignmentMonitor::findPairVertex(const Run3ScoutingMuon& mu1,
                                                                         const Run3ScoutingMuon& mu2,
                                                                         const Run3ScoutingVertexCollection& svs,
                                                                         double& prob) const {
  const Run3ScoutingVertex* best = nullptr;
  prob = -1.;
  const auto& idx2 = mu2.vtxIndx();
  for (const int idx : mu1.vtxIndx()) {
    if (idx < 0 || static_cast<size_t>(idx) >= svs.size())
      continue;
    if (std::find(idx2.begin(), idx2.end(), idx) == idx2.end())
      continue;
    const auto& sv = svs[idx];
    if (!sv.isValidVtx())
      continue;
    const double p = (sv.ndof() > 0) ? TMath::Prob(sv.chi2(), sv.ndof()) : 0.;
    if (p > prob) {
      prob = p;
      best = &sv;
    }
  }
  return best;
}

// ---------------------------------------------------------------------------
// Primary vertex closest (in 3D) to the di-muon vertex
// ---------------------------------------------------------------------------
const Run3ScoutingVertex* ScoutingDiMuonAlignmentMonitor::findClosestVertex(
    const Run3ScoutingVertex& sv, const Run3ScoutingVertexCollection& pvs) const {
  const Run3ScoutingVertex* closest = nullptr;
  double minDist2 = std::numeric_limits<double>::max();
  for (const auto& pv : pvs) {
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
void ScoutingDiMuonAlignmentMonitor::analyze(edm::Event const& iEvent, edm::EventSetup const&) {
  edm::Handle<Run3ScoutingMuonCollection> muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);
  if (!muonHandle.isValid()) {
    edm::LogWarning("ScoutingDiMuonAlignmentMonitor") << "Invalid muon collection handle.";
    return;
  }

  edm::Handle<Run3ScoutingVertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);
  const bool havePV = pvHandle.isValid() && !pvHandle->empty();

  edm::Handle<Run3ScoutingVertexCollection> svHandle;
  iEvent.getByToken(svToken_, svHandle);
  const bool haveSV = svHandle.isValid();

  // ---- muon selection ----
  std::vector<size_t> selected;
  selected.reserve(muonHandle->size());
  for (size_t i = 0; i < muonHandle->size(); ++i) {
    const auto& mu = (*muonHandle)[i];
    if (mu.pt() < minPt_ || std::abs(mu.eta()) > maxEta_)
      continue;
    if (applyMuonID_ && !scoutingDQMUtils::scoutingMuonID(mu))
      continue;
    selected.push_back(i);
  }
  if (selected.size() < 2)
    return;

  for (size_t ii = 0; ii < selected.size(); ++ii) {
    for (size_t jj = ii + 1; jj < selected.size(); ++jj) {
      const auto& mu_a = (*muonHandle)[selected[ii]];
      const auto& mu_b = (*muonHandle)[selected[jj]];
      if (mu_a.charge() * mu_b.charge() >= 0)
        continue;
      const auto& muPlus = mu_a.charge() > 0 ? mu_a : mu_b;
      const auto& muMinus = mu_a.charge() > 0 ? mu_b : mu_a;

      const auto p4 = muonP4(muPlus) + muonP4(muMinus);
      const double mass = p4.mass();
      if (mass < massLimits_.first || mass > massLimits_.second)
        continue;

      // ---- 1. mass bias (all pairs in the window) ----
      hMassVsEtaPlus_->Fill(muPlus.eta(), mass);
      hMassVsEtaMinus_->Fill(muMinus.eta(), mass);
      hMassVsPhiPlus_->Fill(muPlus.phi(), mass);
      hMassVsPhiMinus_->Fill(muMinus.phi(), mass);
      hMassVsDeltaEta_->Fill(muPlus.eta() - muMinus.eta(), mass);
      pMassVsEtaPlus_->Fill(muPlus.eta(), mass);
      pMassVsEtaMinus_->Fill(muMinus.eta(), mass);
      pMassVsPhiPlus_->Fill(muPlus.phi(), mass);
      pMassVsPhiMinus_->Fill(muMinus.phi(), mass);
      pMassVsDeltaEta_->Fill(muPlus.eta() - muMinus.eta(), mass);
      pMassVsEtaPlusEtaMinus_->Fill(muPlus.eta(), muMinus.eta(), mass);

      // ---- vertex-based quantities: need the pair's vertex and a PV ----
      if (!haveSV || !havePV)
        continue;
      double prob = -1.;
      const Run3ScoutingVertex* sv = findPairVertex(muPlus, muMinus, *svHandle, prob);
      if (sv == nullptr || prob < minVtxProb_)
        continue;
      const Run3ScoutingVertex* closestPV = useClosestVertex_ ? findClosestVertex(*sv, *pvHandle) : nullptr;
      const Run3ScoutingVertex& pv = closestPV ? *closestPV : pvHandle->front();

      // ---- 2. impact parameters w.r.t. the PV, by charge ----
      pDxyVsPhiPlus_->Fill(muPlus.phi(), trkDxy(muPlus, pv.x(), pv.y()) * cmToUm);
      pDxyVsPhiMinus_->Fill(muMinus.phi(), trkDxy(muMinus, pv.x(), pv.y()) * cmToUm);
      pDzVsEtaPlus_->Fill(muPlus.eta(), trkDz(muPlus, pv.x(), pv.y(), pv.z()) * cmToUm);
      pDzVsEtaMinus_->Fill(muMinus.eta(), trkDz(muMinus, pv.x(), pv.y(), pv.z()) * cmToUm);

      // ---- 3. PV-independent: both muons w.r.t. their di-muon vertex ----
      const double deltaDxy = (trkDxy(muPlus, sv->x(), sv->y()) - trkDxy(muMinus, sv->x(), sv->y())) * cmToUm;
      const double deltaDz =
          (trkDz(muPlus, sv->x(), sv->y(), sv->z()) - trkDz(muMinus, sv->x(), sv->y(), sv->z())) * cmToUm;
      hDeltaDxy_->Fill(deltaDxy);
      hDeltaDz_->Fill(deltaDz);
      pDeltaDxyVsPhiPlus_->Fill(muPlus.phi(), deltaDxy);
      pDeltaDzVsEtaPlus_->Fill(muPlus.eta(), deltaDz);
      pDeltaDzVsEtaPlusEtaMinus_->Fill(muPlus.eta(), muMinus.eta(), deltaDz);

      // ---- 4. PV -> di-muon vertex displacement ----
      const double dx = sv->x() - pv.x();
      const double dy = sv->y() - pv.y();
      const double dz = sv->z() - pv.z();
      const double pt = p4.pt();
      if (pt <= 0.)
        continue;
      const double lxy = (dx * p4.px() + dy * p4.py()) / pt;     // along the pair p_T
      const double lperp = (-dx * p4.py() + dy * p4.px()) / pt;  // across it
      pVtxDzVsEta_->Fill(p4.eta(), dz * cmToUm);
      pVtxDzVsPhi_->Fill(p4.phi(), dz * cmToUm);
      pLxyVsEta_->Fill(p4.eta(), lxy * cmToUm);
      pLxyVsPhi_->Fill(p4.phi(), lxy * cmToUm);
      pLperpVsPhi_->Fill(p4.phi(), lperp * cmToUm);

      // ---- 5. pseudo-proper decay length and its uncertainty ----
      const double decayLength = lxy * mass / pt;
      const Cov2 c = sumCovXY(*sv, pv);
      const double ux = p4.px() / pt;
      const double uy = p4.py() / pt;
      const double varLxy = ux * ux * c.xx + uy * uy * c.yy + 2. * ux * uy * c.xy;
      const double decayLengthErr = varLxy > 0. ? std::sqrt(varLxy) * mass / pt : 0.;
      hDecayLength_->Fill(decayLength * cmToUm);
      hDecayLengthErr_->Fill(decayLengthErr * cmToUm);
      if (decayLengthErr > 0.)
        hDecayLengthSig_->Fill(decayLength / decayLengthErr);
      if (std::abs(decayLength * cmToUm) < promptWindow_) {
        pDecayLengthVsEta_->Fill(p4.eta(), decayLength * cmToUm);
        pDecayLengthVsPhi_->Fill(p4.phi(), decayLength * cmToUm);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// fillDescriptions
// ---------------------------------------------------------------------------
void ScoutingDiMuonAlignmentMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
  desc.add<double>("promptWindow", 200.)->setComment("um, |l_xy| window for the prompt-peak profiles");

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(ScoutingDiMuonAlignmentMonitor);
