#include "Calibration/HcalCalibAlgos/plugins/HcalPedestalCalibHarvester.h"

#include <map>
#include <memory>
#include <regex>

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "CondFormats/HcalObjects/interface/HcalPedestals.h"
#include "CondFormats/HcalObjects/interface/HcalPedestalWidths.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TH2F.h"
#include "TH1D.h"

namespace {
  HcalSubdetector subdetFromName(std::string const& s) {
    if (s == "HB")
      return HcalBarrel;
    if (s == "HE")
      return HcalEndcap;
    if (s == "HF")
      return HcalForward;
    if (s == "HO")
      return HcalOuter;
    return HcalOther;
  }

  struct ChannelStats {
    HcalDetId did;
    float mean[4] = {0, 0, 0, 0};
    float rms[4] = {0, 0, 0, 0};
    bool valid[4] = {false, false, false, false};
  };
}  // namespace

HcalPedestalCalibHarvester::HcalPedestalCalibHarvester(const edm::ParameterSet& iConfig)
    : minEntriesPerCapId_(iConfig.getUntrackedParameter<double>("minEntriesPerCapId", 10.)),
      dqmDir_(iConfig.getUntrackedParameter<std::string>("dqmDir", "Hcal/PedestalCalib")),
      recordLabel_(iConfig.getUntrackedParameter<std::string>("recordLabel", "")),
      tokenTopology_(esConsumes()) {}

void HcalPedestalCalibHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<double>("minEntriesPerCapId", 10.);
  desc.addUntracked<std::string>("dqmDir", "Hcal/PedestalCalib");
  desc.addUntracked<std::string>("recordLabel", "");
  descriptions.add("hcalPedestalCalibHarvester", desc);
}

void HcalPedestalCalibHarvester::dqmEndRun(DQMStore::IBooker& iBooker,
                                           DQMStore::IGetter& iGetter,
                                           edm::Run const& iRun,
                                           edm::EventSetup const& iSetup) {
  HcalTopology const& topology = iSetup.getData(tokenTopology_);
  HcalPedestals pedestals(&topology, false);
  HcalPedestalWidths widths(&topology, false);

  std::map<uint32_t, ChannelStats> stats;

  static const std::regex nameRe(R"(pedestal_(HB|HE|HF|HO)_ieta(-?\d+)_iphi(\d+)_depth(\d+))");

  for (std::string const sub : {"HB", "HE", "HF", "HO"}) {
    for (MonitorElement* me : iGetter.getAllContents(dqmDir_ + "/" + sub)) {
      std::smatch m;
      std::string const& name = me->getName();
      if (!std::regex_match(name, m, nameRe))
        continue;

      HcalSubdetector subdet = subdetFromName(m[1]);
      int ieta = std::stoi(m[2]);
      int iphi = std::stoi(m[3]);
      int depth = std::stoi(m[4]);
      HcalDetId did(subdet, ieta, iphi, depth);

      TH2F* h2 = me->getTH2F();
      ChannelStats cs;
      cs.did = did;
      for (int capId = 0; capId < 4; ++capId) {
        int xbin = capId + 1;  // bin 1 == capId 0, given axis range (-0.5, 3.5) with 4 bins
        std::unique_ptr<TH1D> proj(h2->ProjectionY((name + "_py").c_str(), xbin, xbin));
        if (proj->GetEntries() >= minEntriesPerCapId_) {
          cs.mean[capId] = proj->GetMean();
          cs.rms[capId] = proj->GetRMS();
          cs.valid[capId] = true;
        }
      }
      stats[did.rawId()] = cs;
    }
  }

  // Interpolate any capId missing enough statistics from eta/phi neighbours (+/-1),
  // same fallback as the original standalone pedestal-extraction script.
  for (auto& [rawId, cs] : stats) {
    for (int capId = 0; capId < 4; ++capId) {
      if (cs.valid[capId])
        continue;
      double sum = 0.;
      int n = 0;
      for (int deta = -1; deta <= 1; ++deta) {
        for (int dphi = -1; dphi <= 1; ++dphi) {
          if (deta == 0 && dphi == 0)
            continue;
          int neta = cs.did.ieta() + deta;
          int nphi = cs.did.iphi() + dphi;
          if (neta == 0)
            neta = (cs.did.ieta() > 0) ? -1 : 1;
          if (nphi == 0)
            nphi = 72;
          if (nphi == 73)
            nphi = 1;
          HcalDetId neighbor(cs.did.subdet(), neta, nphi, cs.did.depth());
          auto it = stats.find(neighbor.rawId());
          if (it != stats.end() && it->second.valid[capId]) {
            sum += it->second.mean[capId];
            ++n;
          }
        }
      }
      if (n > 0) {
        cs.mean[capId] = sum / n;
        edm::LogWarning("HcalPedestalCalibHarvester")
            << "Interpolated capId " << capId << " for channel " << cs.did << " from " << n << " neighbours";
      } else {
        edm::LogWarning("HcalPedestalCalibHarvester")
            << "No statistics and no neighbours for capId " << capId << ", channel " << cs.did;
      }
    }
  }

  for (auto const& [rawId, cs] : stats) {
    HcalPedestal pedItem(rawId, cs.mean[0], cs.mean[1], cs.mean[2], cs.mean[3]);
    pedestals.addValues(pedItem);

    // HcalPedestalWidth stores a full 4x4 capId covariance matrix; we only
    // measured per-capId variance (no cross-capId correlation), so only the
    // diagonal is filled - matches what the original standalone script wrote
    // to DPGfileWidth (zeros off-diagonal).
    HcalPedestalWidth widthItem(rawId);
    for (int capId = 0; capId < 4; ++capId) {
      widthItem.setSigma(capId, capId, cs.rms[capId] * cs.rms[capId]);  // getWidth() returns sqrt(sigma_ii)
    }
    widths.addValues(widthItem);
  }

  edm::Service<cond::service::PoolDBOutputService> dbOutput;
  if (dbOutput.isAvailable()) {
    if (dbOutput->isNewTagRequest("HcalPedestalsRcd"))
      dbOutput->createOneIOV(pedestals, iRun.run(), "HcalPedestalsRcd");
    else
      dbOutput->appendOneIOV(pedestals, iRun.run(), "HcalPedestalsRcd");

    if (dbOutput->isNewTagRequest("HcalPedestalWidthsRcd"))
      dbOutput->createOneIOV(widths, iRun.run(), "HcalPedestalWidthsRcd");
    else
      dbOutput->appendOneIOV(widths, iRun.run(), "HcalPedestalWidthsRcd");
  } else {
    edm::LogError("HcalPedestalCalibHarvester") << "PoolDBOutputService not available - payload not written";
  }
}

DEFINE_FWK_MODULE(HcalPedestalCalibHarvester);
