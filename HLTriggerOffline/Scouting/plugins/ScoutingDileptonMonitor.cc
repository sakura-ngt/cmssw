// system includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// user includes
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "ScoutingDQMUtils.h"

// ---------------------------------------------------------------------------
// Anonymous helpers
// ---------------------------------------------------------------------------
namespace {
  bool checkScoutingID(const Run3ScoutingMuon& obj) { return scoutingDQMUtils::scoutingMuonID(obj); }
  bool checkScoutingID(const Run3ScoutingElectron& obj) { return scoutingDQMUtils::scoutingElectronID(obj); }
}  // namespace

// ---------------------------------------------------------------------------
// scouting helpers - work on plain objects + optional ValueMap
// ---------------------------------------------------------------------------
namespace scouting {

  // Muon: charge is a direct member
  inline int chargeMu(const Run3ScoutingMuon& mu) { return mu.charge(); }

  // Electron: index of the GSF track associated to the electron, -1 if there is none.
  // Use the ValueMap for the best-track index when available, otherwise the first track.
  inline int bestTrackIdxEl(const Run3ScoutingElectron& el,
                            const edm::ValueMap<int>* vmBestIdx,
                            size_t collectionIndex,
                            const edm::Handle<Run3ScoutingElectronCollection>& handle) {
    if (vmBestIdx) {
      edm::Ref<Run3ScoutingElectronCollection> ref(handle, collectionIndex);
      return (*vmBestIdx)[ref];
    }
    return el.trkpt().empty() ? -1 : 0;
  }

  template <typename T>
  math::PtEtaPhiMLorentzVector p4(const T& obj, double mass) {
    return math::PtEtaPhiMLorentzVector(obj.pt(), obj.eta(), obj.phi(), mass);
  }

  // Electron for the J/Psi: Run3ScoutingElectron::eta()/phi() give the supercluster position, far from
  // the direction at the vertex for low-pT electrons. Take the direction from the GSF track (mode
  // estimates when stored, i.e. from the 2024 format on) and keep the SC pT.
  inline math::PtEtaPhiMLorentzVector p4ElTrackDirection(const Run3ScoutingElectron& el, size_t trkIdx) {
    const bool hasMode = !el.trketaMode().empty();
    return math::PtEtaPhiMLorentzVector(el.pt(),
                                        hasMode ? el.trketaMode()[trkIdx] : el.trketa()[trkIdx],
                                        hasMode ? el.trkphiMode()[trkIdx] : el.trkphi()[trkIdx],
                                        scoutingDQMUtils::ELECTRON_MASS);
  }

}  // namespace scouting

// ---------------------------------------------------------------------------
// Monitor class
// ---------------------------------------------------------------------------
class ScoutingDileptonMonitor : public DQMEDAnalyzer {
public:
  explicit ScoutingDileptonMonitor(const edm::ParameterSet&);
  ~ScoutingDileptonMonitor() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(edm::Event const&, edm::EventSetup const&) override;

private:
  // ---- per-collection histogram set -------------------------------------
  struct MassHistos {
    MonitorElement* full{nullptr};
    MonitorElement* zwin{nullptr};
    MonitorElement* jpsiwin{nullptr};
    MonitorElement* barrel{nullptr};
    MonitorElement* endcap{nullptr};
  };

  struct JpsiHistos {
    MonitorElement* lowMass{nullptr};
    MonitorElement* lowMassEBEB{nullptr};
    MonitorElement* jpsiMass{nullptr};
    MonitorElement* jpsiMassEBEB{nullptr};
  };

  // ---- helpers ----------------------------------------------------------
  template <typename T>
  void analyzeCollection(const edm::Event&,
                         const edm::EDGetTokenT<std::vector<T>>&,
                         const edm::EDGetTokenT<edm::ValueMap<int>>*,  // nullptr = no ValueMap (muons)
                         const StringCutObjectSelector<T>&,
                         MassHistos&,
                         bool doEtaSplit,
                         bool useID);

  void fillPairs(const std::vector<size_t>& indices,
                 const edm::Handle<std::vector<Run3ScoutingMuon>>& handle,
                 const edm::ValueMap<int>* /*unused*/,
                 MassHistos&,
                 bool doEtaSplit);

  void fillPairs(const std::vector<size_t>& indices,
                 const edm::Handle<std::vector<Run3ScoutingElectron>>& handle,
                 const edm::ValueMap<int>* vmBestIdx,
                 MassHistos&,
                 bool doEtaSplit);

  bool passJpsiElectronTrigger(const edm::Event&) const;
  void analyzeJpsiElectrons(const edm::Event&);

  // ---- configuration ----------------------------------------------------
  const std::string outputInternalPath_;
  const double massMin_, massMax_;
  const int massBins_;
  const double zMin_, zMax_;
  const double jpsiMin_, jpsiMax_;
  const double barrelEta_;

  // ---- muons (Vtx) ------------------------------------------------------
  const bool doMuons_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muonToken_;
  const StringCutObjectSelector<Run3ScoutingMuon> muonCut_;
  MassHistos muonHistos_;
  const bool muonID_;

  // ---- muons (No Vtx) ---------------------------------------------------
  const bool doMuonsNoVtx_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muonNoVtxToken_;
  MassHistos muonNoVtxHistos_;

  // ---- electrons --------------------------------------------------------
  const bool doElectrons_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron>> electronToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> vmBestTrackIndexToken_;
  const StringCutObjectSelector<Run3ScoutingElectron> electronCut_;
  MassHistos electronHistos_;
  const bool electronID_;

  // ---- J/Psi -> ee (selection of the EGM scouting studies) ----------------
  const bool doJpsiElectrons_;
  const edm::EDGetTokenT<double> rhoToken_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  const std::vector<std::string> jpsiEleTriggers_;
  const double jpsiElePtMin_, jpsiEleTrackEtaMax_;
  const scoutingDQMUtils::ElectronWP jpsiEleWP_;
  const bool jpsiEleUseIso_;
  const int jpsiEleLowMassBins_;
  const double jpsiEleLowMassMax_;
  const int jpsiEleMassBins_;
  const double jpsiEleMassMin_, jpsiEleMassMax_;
  JpsiHistos jpsiElectronHistos_;
};

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
ScoutingDileptonMonitor::ScoutingDileptonMonitor(const edm::ParameterSet& iConfig)
    : outputInternalPath_{iConfig.getParameter<std::string>("OutputInternalPath")},
      massMin_(iConfig.getParameter<double>("massMin")),
      massMax_(iConfig.getParameter<double>("massMax")),
      massBins_(iConfig.getParameter<int>("massBins")),
      zMin_(iConfig.getParameter<double>("zMassMin")),
      zMax_(iConfig.getParameter<double>("zMassMax")),
      jpsiMin_(iConfig.getParameter<double>("jpsiMassMin")),
      jpsiMax_(iConfig.getParameter<double>("jpsiMassMax")),
      barrelEta_(iConfig.getParameter<double>("barrelEta")),
      doMuons_(iConfig.getParameter<bool>("doMuons")),
      muonToken_(consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      muonCut_(iConfig.getParameter<std::string>("muonCut")),
      muonID_(iConfig.getParameter<bool>("muonID")),
      doMuonsNoVtx_(iConfig.getParameter<bool>("doMuonsNoVtx")),
      muonNoVtxToken_(consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muonsNoVtx"))),
      doElectrons_(iConfig.getParameter<bool>("doElectrons")),
      electronToken_(consumes<std::vector<Run3ScoutingElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
      vmBestTrackIndexToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("vmBestTrackIndex"))),
      electronCut_(iConfig.getParameter<std::string>("electronCut")),
      electronID_(iConfig.getParameter<bool>("electronID")),
      doJpsiElectrons_(iConfig.getParameter<bool>("doJpsiElectrons")),
      rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
      triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      jpsiEleTriggers_(iConfig.getParameter<std::vector<std::string>>("jpsiEleTriggers")),
      jpsiElePtMin_(iConfig.getParameter<double>("jpsiElePtMin")),
      jpsiEleTrackEtaMax_(iConfig.getParameter<double>("jpsiEleTrackEtaMax")),
      jpsiEleWP_(scoutingDQMUtils::electronWPFromString(iConfig.getParameter<std::string>("jpsiEleWP"))),
      jpsiEleUseIso_(iConfig.getParameter<bool>("jpsiEleUseIso")),
      jpsiEleLowMassBins_(iConfig.getParameter<int>("jpsiEleLowMassBins")),
      jpsiEleLowMassMax_(iConfig.getParameter<double>("jpsiEleLowMassMax")),
      jpsiEleMassBins_(iConfig.getParameter<int>("jpsiEleMassBins")),
      jpsiEleMassMin_(iConfig.getParameter<double>("jpsiEleMassMin")),
      jpsiEleMassMax_(iConfig.getParameter<double>("jpsiEleMassMax")) {}

// ---------------------------------------------------------------------------
// bookHistograms
// ---------------------------------------------------------------------------
void ScoutingDileptonMonitor::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const&, edm::EventSetup const&) {
  ibooker.setCurrentFolder(outputInternalPath_);

  auto bookSet = [&](const std::string& name, MassHistos& h, bool splitEta) {
    h.full = ibooker.book1D(
        name + "_mass", name + " opposite-charge invariant mass;M [GeV];Events", massBins_, massMin_, massMax_);
    h.zwin = ibooker.book1D(name + "_zMass", name + " Z window;M [GeV];Events", massBins_, zMin_, zMax_);
    h.jpsiwin =
        ibooker.book1D(name + "_jpsiMass", name + " J/Psi window;M [GeV];Events", massBins_, jpsiMin_, jpsiMax_);
    if (splitEta) {
      h.barrel = ibooker.book1D(name + "_barrelMass", name + " barrel;M [GeV];Events", massBins_, massMin_, massMax_);
      h.endcap = ibooker.book1D(name + "_endcapMass", name + " endcap;M [GeV];Events", massBins_, massMin_, massMax_);
    }
  };

  if (doMuons_)
    bookSet("muons", muonHistos_, true);
  if (doMuonsNoVtx_)
    bookSet("muonsNoVtx", muonNoVtxHistos_, true);
  if (doElectrons_)
    bookSet("electrons", electronHistos_, true);

  if (doElectrons_ && doJpsiElectrons_) {
    const std::string name = "electronsJpsi";
    auto& h = jpsiElectronHistos_;
    // ROOT TLatex titles; y axis carries the bin width
    auto bookMass = [&](const std::string& suffix, const std::string& title, int bins, double min, double max) {
      const std::string axes =
          std::format(";M_{{e^{{+}}e^{{-}}}} [GeV];Events / {:g} GeV", std::round((max - min) / bins * 1e4) / 1e4);
      return ibooker.book1D(name + suffix, title + axes, bins, min, max);
    };
    const std::string pairTitle = "Scouting e^{+}e^{-}, leading opposite-charge pair";
    const std::string jpsiTitle = "Scouting J/#psi #rightarrow e^{+}e^{-}";
    const std::string ebebTitle = " (both e in ECAL barrel)";
    h.lowMass = bookMass("_lowMass", pairTitle, jpsiEleLowMassBins_, 0., jpsiEleLowMassMax_);
    h.lowMassEBEB = bookMass("_lowMassEBEB", pairTitle + ebebTitle, jpsiEleLowMassBins_, 0., jpsiEleLowMassMax_);
    h.jpsiMass = bookMass("_jpsiMass", jpsiTitle, jpsiEleMassBins_, jpsiEleMassMin_, jpsiEleMassMax_);
    h.jpsiMassEBEB =
        bookMass("_jpsiMassEBEB", jpsiTitle + ebebTitle, jpsiEleMassBins_, jpsiEleMassMin_, jpsiEleMassMax_);
  }
}

// ---------------------------------------------------------------------------
// analyze
// ---------------------------------------------------------------------------
void ScoutingDileptonMonitor::analyze(edm::Event const& iEvent, edm::EventSetup const&) {
  if (doMuons_) {
    edm::LogInfo("ScoutingDileptonMonitor") << "doing muons";
    analyzeCollection(iEvent, muonToken_, nullptr, muonCut_, muonHistos_, true, muonID_);
  }
  if (doMuonsNoVtx_) {
    edm::LogInfo("ScoutingDileptonMonitor") << "doing muons NoVtx";
    analyzeCollection(iEvent, muonNoVtxToken_, nullptr, muonCut_, muonNoVtxHistos_, true, muonID_);
  }
  if (doElectrons_) {
    edm::LogInfo("ScoutingDileptonMonitor") << "doing electrons";
    analyzeCollection(
        iEvent, electronToken_, &vmBestTrackIndexToken_, electronCut_, electronHistos_, true, electronID_);
    if (doJpsiElectrons_)
      analyzeJpsiElectrons(iEvent);
  }
}

// ---------------------------------------------------------------------------
// analyzeCollection
// ---------------------------------------------------------------------------
template <typename T>
void ScoutingDileptonMonitor::analyzeCollection(const edm::Event& iEvent,
                                                const edm::EDGetTokenT<std::vector<T>>& token,
                                                const edm::EDGetTokenT<edm::ValueMap<int>>* bestTrackIndexToken,
                                                const StringCutObjectSelector<T>& cut,
                                                MassHistos& histos,
                                                bool doEtaSplit,
                                                bool useID) {
  edm::Handle<std::vector<T>> handle;
  iEvent.getByToken(token, handle);

  if (!handle.isValid()) {
    edm::LogWarning("ScoutingDileptonMonitor") << "Invalid Handle!";
    return;
  }

  const edm::ValueMap<int>* vmBestIdx = nullptr;
  if (bestTrackIndexToken != nullptr)
    vmBestIdx = &iEvent.get(*bestTrackIndexToken);

  edm::LogInfo("ScoutingDileptonMonitor") << "collection size: " << handle->size();

  std::vector<size_t> selected;
  selected.reserve(handle->size());
  for (size_t i = 0; i < handle->size(); ++i) {
    const auto& obj = (*handle)[i];
    if (!cut(obj))
      continue;
    if (useID && !checkScoutingID(obj))
      continue;
    selected.push_back(i);
  }

  fillPairs(selected, handle, vmBestIdx, histos, doEtaSplit);
}

// ---------------------------------------------------------------------------
// fillPairs - muon overload
// ---------------------------------------------------------------------------
void ScoutingDileptonMonitor::fillPairs(const std::vector<size_t>& indices,
                                        const edm::Handle<std::vector<Run3ScoutingMuon>>& handle,
                                        const edm::ValueMap<int>* /*unused*/,
                                        MassHistos& histos,
                                        bool doEtaSplit) {
  const size_t n = indices.size();
  edm::LogInfo("ScoutingDileptonMonitor") << "muon lepton size: " << n;

  for (size_t ii = 0; ii < n; ++ii) {
    for (size_t jj = ii + 1; jj < n; ++jj) {
      const auto& mu_i = (*handle)[indices[ii]];
      const auto& mu_j = (*handle)[indices[jj]];

      if (scouting::chargeMu(mu_i) * scouting::chargeMu(mu_j) >= 0)
        continue;

      const auto p4 = scouting::p4(mu_i, scoutingDQMUtils::MUON_MASS) + scouting::p4(mu_j, scoutingDQMUtils::MUON_MASS);
      const double mass = p4.mass();

      histos.full->Fill(mass);
      if (mass > zMin_ && mass < zMax_)
        histos.zwin->Fill(mass);
      if (mass > jpsiMin_ && mass < jpsiMax_)
        histos.jpsiwin->Fill(mass);
      if (doEtaSplit) {
        const bool barrel = std::abs(mu_i.eta()) < barrelEta_ && std::abs(mu_j.eta()) < barrelEta_;
        (barrel ? histos.barrel : histos.endcap)->Fill(mass);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// fillPairs - electron overload
// ---------------------------------------------------------------------------
void ScoutingDileptonMonitor::fillPairs(const std::vector<size_t>& indices,
                                        const edm::Handle<std::vector<Run3ScoutingElectron>>& handle,
                                        const edm::ValueMap<int>* vmBestIdx,
                                        MassHistos& histos,
                                        bool doEtaSplit) {
  // keep only electrons with an associated track: needed for the charge
  std::vector<std::pair<size_t, int>> withTrack;  // (collection index, best track index)
  withTrack.reserve(indices.size());
  for (const size_t idx : indices) {
    const int trkIdx = scouting::bestTrackIdxEl((*handle)[idx], vmBestIdx, idx, handle);
    if (trkIdx >= 0)
      withTrack.emplace_back(idx, trkIdx);
  }

  const size_t n = withTrack.size();
  edm::LogInfo("ScoutingDileptonMonitor") << "electron lepton size: " << n;

  for (size_t ii = 0; ii < n; ++ii) {
    for (size_t jj = ii + 1; jj < n; ++jj) {
      const auto [idx_i, trk_i] = withTrack[ii];
      const auto [idx_j, trk_j] = withTrack[jj];
      const auto& el_i = (*handle)[idx_i];
      const auto& el_j = (*handle)[idx_j];

      if (el_i.trkcharge()[trk_i] * el_j.trkcharge()[trk_j] >= 0)
        continue;

      const auto p4 =
          scouting::p4(el_i, scoutingDQMUtils::ELECTRON_MASS) + scouting::p4(el_j, scoutingDQMUtils::ELECTRON_MASS);
      const double mass = p4.mass();

      histos.full->Fill(mass);
      if (mass > zMin_ && mass < zMax_)
        histos.zwin->Fill(mass);
      if (mass > jpsiMin_ && mass < jpsiMax_)
        histos.jpsiwin->Fill(mass);
      if (doEtaSplit) {
        const bool barrel = std::abs(el_i.eta()) < barrelEta_ && std::abs(el_j.eta()) < barrelEta_;
        (barrel ? histos.barrel : histos.endcap)->Fill(mass);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// passJpsiElectronTrigger - any path starting with one of jpsiEleTriggers fired
// ---------------------------------------------------------------------------
bool ScoutingDileptonMonitor::passJpsiElectronTrigger(const edm::Event& iEvent) const {
  // if (jpsiEleTriggers_.empty())
  //   return true;

  // const auto& triggerResults = iEvent.getHandle(triggerResultsToken_);
  // if (!triggerResults.isValid()) {
  //   edm::LogWarning("ScoutingDileptonMonitor") << "Invalid TriggerResults handle, skipping J/Psi -> ee";
  //   return false;
  // }

  // const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
  // for (unsigned int i = 0; i < triggerNames.size(); ++i) {
  //   if (!triggerResults->accept(i))
  //     continue;
  //   for (const auto& pattern : jpsiEleTriggers_) {
  //     if (triggerNames.triggerName(i).find(pattern) == 0)
  //       return true;
  //   }
  // }
  return true;
}

// ---------------------------------------------------------------------------
// analyzeJpsiElectrons - mirrors the dielectron mass of the EGM scouting studies:
// pT(SC) > jpsiElePtMin, best GSF track (scoutingElectronGsfTrackIdx) with
// |eta| < jpsiEleTrackEtaMax, cut-based ID, two leading electrons with opposite charge,
// direction from the GSF track
// ---------------------------------------------------------------------------
void ScoutingDileptonMonitor::analyzeJpsiElectrons(const edm::Event& iEvent) {
  const auto& handle = iEvent.getHandle(electronToken_);
  if (!handle.isValid() || !passJpsiElectronTrigger(iEvent))
    return;

  const auto& rhoHandle = iEvent.getHandle(rhoToken_);
  const double rho = rhoHandle.isValid() ? *rhoHandle : 0.;

  std::vector<std::pair<const Run3ScoutingElectron*, size_t>> selected;  // (electron, best track index)
  for (const auto& el : *handle) {
    if (el.pt() <= jpsiElePtMin_)
      continue;
    size_t trkIdx = 0;
    if (!scoutingDQMUtils::scoutingElectronGsfTrackIdx(el, trkIdx))
      continue;
    if (std::abs(el.trketa()[trkIdx]) >= jpsiEleTrackEtaMax_)
      continue;
    if (!scoutingDQMUtils::scoutingElectronCutBasedID(el, rho, jpsiEleWP_, jpsiEleUseIso_))
      continue;
    selected.emplace_back(&el, trkIdx);
  }

  edm::LogInfo("ScoutingDileptonMonitor") << "J/Psi electron size: " << selected.size();
  if (selected.size() < 2)
    return;

  std::partial_sort(selected.begin(), selected.begin() + 2, selected.end(), [](const auto& a, const auto& b) {
    return a.first->pt() > b.first->pt();
  });
  const auto& [el_1, trk_1] = selected[0];
  const auto& [el_2, trk_2] = selected[1];

  if (el_1->trkcharge()[trk_1] * el_2->trkcharge()[trk_2] >= 0)
    return;

  const double mass = (scouting::p4ElTrackDirection(*el_1, trk_1) + scouting::p4ElTrackDirection(*el_2, trk_2)).mass();
  const bool ebeb = std::abs(el_1->eta()) < barrelEta_ && std::abs(el_2->eta()) < barrelEta_;

  auto& h = jpsiElectronHistos_;
  h.lowMass->Fill(mass);
  h.jpsiMass->Fill(mass);
  if (ebeb) {
    h.lowMassEBEB->Fill(mass);
    h.jpsiMassEBEB->Fill(mass);
  }
}

// ---------------------------------------------------------------------------
// fillDescriptions
// ---------------------------------------------------------------------------
void ScoutingDileptonMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("OutputInternalPath", "HLT/ScoutingOffline/DiLepton");
  desc.add<edm::InputTag>("muons", edm::InputTag("hltScoutingMuonPackerVtx"));
  desc.add<edm::InputTag>("muonsNoVtx", edm::InputTag("hltScoutingMuonPackerNoVtx"));
  desc.add<edm::InputTag>("electrons", edm::InputTag("hltScoutingEgammaPacker"));
  desc.add<edm::InputTag>("vmBestTrackIndex",
                          edm::InputTag("run3ScoutingElectronBestTrack", "Run3ScoutingElectronBestTrackIndex"));
  desc.add<bool>("doMuons", true);
  desc.add<bool>("doMuonsNoVtx", true);
  desc.add<bool>("doElectrons", true);
  desc.add<std::string>("muonCut", "");
  desc.add<std::string>("electronCut", "");
  desc.add<int>("massBins", 120);
  desc.add<double>("massMin", 0.0);
  desc.add<double>("massMax", 200.0);
  desc.add<double>("zMassMin", 70.0);
  desc.add<double>("zMassMax", 110.0);
  desc.add<double>("jpsiMassMin", 2.6);
  desc.add<double>("jpsiMassMax", 3.5);
  desc.add<double>("barrelEta", 1.479);
  desc.add<bool>("muonID", true);
  desc.add<bool>("electronID", true);
  desc.add<bool>("doJpsiElectrons", true)
      ->setComment("J/Psi -> ee with the selection of the EGM scouting studies (needs doElectrons)");
  desc.add<edm::InputTag>("rho", edm::InputTag("hltScoutingPFPacker", "rho"));
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults", "", "HLT"));
  desc.add<std::vector<std::string>>("jpsiEleTriggers", {"DST_PFScouting_DoubleEG"})
      ->setComment("path-name prefixes, at least one must have fired; empty = no requirement");
  desc.add<double>("jpsiElePtMin", 5.0);
  desc.add<double>("jpsiEleTrackEtaMax", 2.5);
  desc.add<std::string>("jpsiEleWP", "Loose")->setComment("cut-based ID: Veto, Loose, Medium or Tight");
  desc.add<bool>("jpsiEleUseIso", false);
  desc.add<int>("jpsiEleLowMassBins", 200);
  desc.add<double>("jpsiEleLowMassMax", 20.0);
  desc.add<int>("jpsiEleMassBins", 100);
  desc.add<double>("jpsiEleMassMin", 2.0);
  desc.add<double>("jpsiEleMassMax", 4.0);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(ScoutingDileptonMonitor);
