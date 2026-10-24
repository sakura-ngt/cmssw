#ifndef HLTriggerOffline_Scouting_ScoutingDQMUtils_h
#define HLTriggerOffline_Scouting_ScoutingDQMUtils_h

#include <algorithm>
#include <cmath>
#include <string>

#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace scoutingDQMUtils {

  // Constants
  static constexpr double MUON_MASS = 0.105658;      // Muon mass in GeV
  static constexpr double ELECTRON_MASS = 0.000511;  // Electron mass in GeV
  static constexpr double ELE_etaEB = 1.479;         // Eta restriction to barrel (for electrons)

  // trigs
  inline double computePtFromEnergyMassEta(double energy, double mass, double eta) {
    double transverseEnergy = std::sqrt(energy * energy - mass * mass);
    return transverseEnergy / std::cosh(eta);
  }

  // scouting muons ID
  inline const bool scoutingMuonID(const Run3ScoutingMuon& mu) {
    if (mu.pt() < 3.0)
      return false;
    if (std::abs(mu.eta()) > 2.4)
      return false;
    if (mu.normalizedChi2() > 3)
      return false;
    return true;
  }

  // scouting electron IB
  inline const bool scoutingElectronID(const Run3ScoutingElectron& el) {
    bool isEB = (std::abs(el.eta()) < scoutingDQMUtils::ELE_etaEB);
    if (isEB) {
      if (el.sigmaIetaIeta() > 0.0146)
        return false;
      if (el.hOverE() > 0.2)
        return false;
      if (std::abs(el.ooEMOop()) > 0.0641)
        return false;
      if (std::abs(el.dEtaIn()) > 0.0087)
        return false;
      if (std::abs(el.dPhiIn()) > 0.2653)
        return false;
      return true;

    } else {
      if (el.sigmaIetaIeta() > 0.0305)
        return false;
      if (el.hOverE() > 0.2)
        return false;
      if (std::abs(el.dEtaIn()) > 1.5994)
        return false;
      if (std::abs(el.dPhiIn()) > 0.2728)
        return false;
      if (std::abs(el.ooEMOop()) > 0.0681)
        return false;
      return true;
    }
  }

  inline bool scoutingElectronGsfTrackID(const Run3ScoutingElectron& el, size_t trackIdx) {
    if (trackIdx > el.trkpt().size())
      edm::LogError("ScoutingDQMUtils") << "Invalid track index for electron: Exceeds the number of tracks";

    math::PtEtaPhiMLorentzVector particleSC(el.pt(), el.eta(), el.phi(), scoutingDQMUtils::ELECTRON_MASS);
    math::PtEtaPhiMLorentzVector particleTrk(
        el.trkpt()[trackIdx], el.trketa()[trackIdx], el.trkphi()[trackIdx], scoutingDQMUtils::ELECTRON_MASS);

    double scEnergy = particleSC.energy();
    double trkEnergy = particleTrk.energy();
    double relEnergyDiff = std::abs(scEnergy - trkEnergy) / scEnergy;
    double dPhi = deltaPhi(particleSC.phi(), particleTrk.phi());

    bool isEB = (std::abs(el.eta()) < scoutingDQMUtils::ELE_etaEB);
    if (isEB) {
      if (el.trkpt()[trackIdx] < 12)
        return false;
      if (relEnergyDiff > 1)
        return false;
      if (dPhi > 0.06)
        return false;
      if (el.trkchi2overndf()[trackIdx] > 3)
        return false;
      return true;
    } else {
      if (el.trkpt()[trackIdx] < 12)
        return false;
      if (relEnergyDiff > 1)
        return false;
      if (dPhi > 0.06)
        return false;
      if (el.trkchi2overndf()[trackIdx] > 2)
        return false;
      return true;
    }
  }

  inline bool scoutingElectronGsfTrackIdx(const Run3ScoutingElectron& el, size_t& trackIdx) {
    bool foundGoodGsfTrkIdx = false;
    for (size_t i = 0; i < el.trkpt().size(); ++i) {
      if (scoutingDQMUtils::scoutingElectronGsfTrackID(el, i)) {
        if (!foundGoodGsfTrkIdx) {
          foundGoodGsfTrkIdx = true;
          trackIdx = i;
        } else {
          double relPtDiff = fabs(el.trkpt()[i] - el.pt()) / el.pt();
          double relPtDiffOld = fabs(el.trkpt()[trackIdx] - el.pt()) / el.pt();
          if (relPtDiff < relPtDiffOld)
            trackIdx = i;
        }
      }
    }
    return foundGoodGsfTrkIdx;
  }

  // Cut-based scouting electron ID (122X-tuned working points) of the EGM scouting studies
  enum class ElectronWP { Veto = 0, Loose, Medium, Tight };

  inline ElectronWP electronWPFromString(const std::string& wp) {
    if (wp == "Veto")
      return ElectronWP::Veto;
    if (wp == "Loose")
      return ElectronWP::Loose;
    if (wp == "Medium")
      return ElectronWP::Medium;
    if (wp == "Tight")
      return ElectronWP::Tight;
    throw cms::Exception("Configuration")
        << "Unknown electron working point '" << wp << "', allowed: Veto, Loose, Medium, Tight";
  }

  struct ElectronCutBasedCuts {
    float sieie;
    float dEtaIn;
    float dPhiIn;
    float hoeC0, hoeCE, hoeCr;  // H/E < C0 + CE/E_SC + Cr*rho/E_SC
    float isoC0, isoCpt;        // relIso < C0 + Cpt/pT
    float ooEmooP;
    int missingHits;
  };

  // [working point][0 = EB, 1 = EE]
  inline constexpr ElectronCutBasedCuts electronCutBasedCuts[4][2] = {
      {{0.0117, 0.0071, 0.208, 0.05, 1.28, 0.0422, 0.406, 0.535, 0.178, 2},
       {0.0298, 0.0173, 0.234, 0.05, 2.3, 0.262, 0.342, 0.519, 0.137, 3}},
      {{0.0107, 0.00691, 0.175, 0.05, 1.28, 0.0422, 0.194, 0.535, 0.138, 1},
       {0.0275, 0.0121, 0.228, 0.05, 2.3, 0.262, 0.184, 0.519, 0.127, 1}},
      {{0.0103, 0.00481, 0.127, 0.0241, 1.28, 0.0422, 0.0837, 0.535, 0.0966, 1},
       {0.0272, 0.00951, 0.221, 0.05, 2.3, 0.262, 0.0741, 0.519, 0.111, 1}},
      {{0.0101, 0.00411, 0.116, 0.02, 1.16, 0.0422, 0.0388, 0.535, 0.023, 1},
       {0.027, 0.00938, 0.164, 0.02, 0.5, 0.262, 0.0544, 0.519, 0.018, 1}}};

  inline double electronEffectiveArea(double absEta) {
    if (absEta < 1.0)
      return 0.1243;
    if (absEta < 1.479)
      return 0.1458;
    if (absEta < 2.0)
      return 0.0992;
    if (absEta < 2.2)
      return 0.0794;
    if (absEta < 2.3)
      return 0.0762;
    if (absEta < 2.4)
      return 0.0766;
    return 0.1003;
  }

  inline bool scoutingElectronCutBasedID(const Run3ScoutingElectron& el, double rho, ElectronWP wp, bool useIso) {
    const double absEta = std::abs(el.eta());
    const auto& c = electronCutBasedCuts[static_cast<int>(wp)][absEta <= scoutingDQMUtils::ELE_etaEB ? 0 : 1];
    const double eSC = el.rawEnergy();
    if (el.sigmaIetaIeta() >= c.sieie)
      return false;
    if (std::abs(el.dEtaIn()) >= c.dEtaIn)
      return false;
    if (std::abs(el.dPhiIn()) >= c.dPhiIn)
      return false;
    if (el.hOverE() >= c.hoeC0 + c.hoeCE / eSC + c.hoeCr * rho / eSC)
      return false;
    if (std::abs(el.ooEMOop()) >= c.ooEmooP)
      return false;
    if (el.missingHits() > c.missingHits)
      return false;
    if (useIso) {
      const double neutralIso = el.ecalIso() + el.hcalIso() - rho * electronEffectiveArea(absEta);
      const double relIso = (std::max(neutralIso, 0.) + el.trackIso()) / el.pt();
      if (relIso >= c.isoC0 + c.isoCpt / el.pt())
        return false;
    }
    return true;
  }

  inline bool hasPatternInHLTPath(const edm::TriggerNames& triggerNames, const std::string& pattern) {
    for (unsigned int i = 0; i < triggerNames.size(); ++i) {
      const std::string& triggerName = triggerNames.triggerName(i);

      // Check if triggerName starts with the specified prefix
      if (triggerName.find(pattern) == 0) {  // Position 0 means it starts with 'prefix'
        return true;                         // Pattern match found
      }
    }
    return false;  // No match found
  }
}  // namespace scoutingDQMUtils

#endif  //
