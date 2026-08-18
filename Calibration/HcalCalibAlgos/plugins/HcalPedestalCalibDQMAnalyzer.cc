#include "Calibration/HcalCalibAlgos/plugins/HcalPedestalCalibDQMAnalyzer.h"

#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"  // VERIFY: only if needed; adc2fC path below may not require it
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/MakerMacros.h"

namespace {
  // The three digi types disagree on whether the DetId accessor is detid() or id() —
  // matches HcalDigiSortedTableProducer's own per-subdet handling.
  HcalDetId digiDetId(QIE11DataFrame const& digi) { return digi.detid(); }
  HcalDetId digiDetId(QIE10DataFrame const& digi) { return digi.detid(); }
  HcalDetId digiDetId(HODataFrame const& digi) { return digi.id(); }

  std::string subdetName(HcalSubdetector subdet) {
    switch (subdet) {
      case HcalBarrel:
        return "HB";
      case HcalEndcap:
        return "HE";
      case HcalForward:
        return "HF";
      case HcalOuter:
        return "HO";
      default:
        return "Unknown";
    }
  }
}  // namespace

HcalPedestalCalibDQMAnalyzer::HcalPedestalCalibDQMAnalyzer(const edm::ParameterSet& iConfig)
    : tokenHBDetIdList_(consumes<edm::InRun>(
          iConfig.getUntrackedParameter<edm::InputTag>("HBDetIdList", edm::InputTag("hcalDetIdTable", "HBDetIdList")))),
      tokenHEDetIdList_(consumes<edm::InRun>(
          iConfig.getUntrackedParameter<edm::InputTag>("HEDetIdList", edm::InputTag("hcalDetIdTable", "HEDetIdList")))),
      tokenHFDetIdList_(consumes<edm::InRun>(
          iConfig.getUntrackedParameter<edm::InputTag>("HFDetIdList", edm::InputTag("hcalDetIdTable", "HFDetIdList")))),
      tokenHODetIdList_(consumes<edm::InRun>(
          iConfig.getUntrackedParameter<edm::InputTag>("HODetIdList", edm::InputTag("hcalDetIdTable", "HODetIdList")))),
      tokenQIE11_(consumes<QIE11DigiCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>("tagQIE11", edm::InputTag("hcalDigis")))),
      tokenQIE10_(consumes<QIE10DigiCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>("tagQIE10", edm::InputTag("hcalDigis")))),
      tokenHO_(consumes<HODigiCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>("tagHO", edm::InputTag("hcalDigis")))),
      tokenDbService_(esConsumes()),
      fcAxisBins_(iConfig.getUntrackedParameter<int>("fcAxisBins", 2000)),
      fcAxisMax_(iConfig.getUntrackedParameter<double>("fcAxisMax", 2000.0)) {}

void HcalPedestalCalibDQMAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("HBDetIdList", edm::InputTag("hcalDetIdTable", "HBDetIdList"));
  desc.addUntracked<edm::InputTag>("HEDetIdList", edm::InputTag("hcalDetIdTable", "HEDetIdList"));
  desc.addUntracked<edm::InputTag>("HFDetIdList", edm::InputTag("hcalDetIdTable", "HFDetIdList"));
  desc.addUntracked<edm::InputTag>("HODetIdList", edm::InputTag("hcalDetIdTable", "HODetIdList"));
  desc.addUntracked<edm::InputTag>("tagQIE11", edm::InputTag("hcalDigis"));
  desc.addUntracked<edm::InputTag>("tagQIE10", edm::InputTag("hcalDigis"));
  desc.addUntracked<edm::InputTag>("tagHO", edm::InputTag("hcalDigis"));
  desc.addUntracked<int>("fcAxisBins", 2000);
  desc.addUntracked<double>("fcAxisMax", 2000.0);
  descriptions.add("hcalPedestalCalibDQMAnalyzer", desc);
}

void HcalPedestalCalibDQMAnalyzer::bookSubdet(DQMStore::IBooker& iBooker,
                                              edm::Run const& iRun,
                                              edm::EDGetTokenT<std::vector<HcalDetId>> const& token,
                                              HcalSubdetector subdet,
                                              std::string const& folderSuffix) {
  edm::Handle<std::vector<HcalDetId>> detIds;
  iRun.getByToken(token, detIds);
  if (!detIds.isValid())
    return;

  iBooker.setCurrentFolder("Hcal/PedestalCalib/" + folderSuffix);
  for (HcalDetId const& did : *detIds) {
    // Name encodes (subdet, ieta, iphi, depth) so the harvester can recover the
    // channel identity from the booked-histogram name alone, with no ES access needed.
    std::string name = "pedestal_" + subdetName(subdet) + "_ieta" + std::to_string(did.ieta()) + "_iphi" +
                       std::to_string(did.iphi()) + "_depth" + std::to_string(did.depth());
    MonitorElement* me =
        iBooker.book2D(name, "Pedestal charge vs capId;capId;Charge [fC]", 4, -0.5, 3.5, fcAxisBins_, 0., fcAxisMax_);
    pedestalByCapId_[did.rawId()] = me;
  }
}

void HcalPedestalCalibDQMAnalyzer::bookHistograms(DQMStore::IBooker& iBooker,
                                                  edm::Run const& iRun,
                                                  edm::EventSetup const&) {
  pedestalByCapId_.clear();
  bookSubdet(iBooker, iRun, tokenHBDetIdList_, HcalBarrel, "HB");
  bookSubdet(iBooker, iRun, tokenHEDetIdList_, HcalEndcap, "HE");
  bookSubdet(iBooker, iRun, tokenHFDetIdList_, HcalForward, "HF");
  bookSubdet(iBooker, iRun, tokenHODetIdList_, HcalOuter, "HO");
}

void HcalPedestalCalibDQMAnalyzer::fillHBHE(QIE11DigiCollection const& digis,
                                            HcalSubdetector subdet,
                                            HcalDbService const& dbService) {
  for (auto const& raw : digis) {
    const QIE11DataFrame digi = static_cast<const QIE11DataFrame>(raw);
    HcalDetId const did = digiDetId(digi);
    if (did.subdet() != subdet)
      continue;

    auto it = pedestalByCapId_.find(did.rawId());
    if (it == pedestalByCapId_.end())
      continue;

    const HcalQIECoder* channelCoder = dbService.getHcalCoder(did);
    const HcalQIEShape* shape = dbService.getHcalShape(channelCoder);
    HcalCoderDb coder(*channelCoder, *shape);
    CaloSamples samples;
    coder.adc2fC(digi, samples);

    for (int ts = 0; ts < samples.size(); ++ts) {
      it->second->Fill(digi[ts].capid(), samples[ts]);
    }
  }
}

void HcalPedestalCalibDQMAnalyzer::fillHF(QIE10DigiCollection const& digis, HcalDbService const& dbService) {
  for (auto const& raw : digis) {
    const QIE10DataFrame digi = static_cast<const QIE10DataFrame>(raw);
    HcalDetId const did = digiDetId(digi);
    if (did.subdet() != HcalForward)
      continue;

    auto it = pedestalByCapId_.find(did.rawId());
    if (it == pedestalByCapId_.end())
      continue;

    const HcalQIECoder* channelCoder = dbService.getHcalCoder(did);
    const HcalQIEShape* shape = dbService.getHcalShape(channelCoder);
    HcalCoderDb coder(*channelCoder, *shape);
    CaloSamples samples;
    coder.adc2fC(digi, samples);

    for (int ts = 0; ts < samples.size(); ++ts) {
      it->second->Fill(digi[ts].capid(), samples[ts]);
    }
  }
}

void HcalPedestalCalibDQMAnalyzer::fillHO(HODigiCollection const& digis, HcalDbService const& dbService) {
  for (auto const& raw : digis) {
    const HODataFrame digi = static_cast<const HODataFrame>(raw);
    HcalDetId const did = digiDetId(digi);
    if (did.subdet() != HcalOuter)
      continue;

    auto it = pedestalByCapId_.find(did.rawId());
    if (it == pedestalByCapId_.end())
      continue;

    const HcalQIECoder* channelCoder = dbService.getHcalCoder(did);
    const HcalQIEShape* shape = dbService.getHcalShape(channelCoder);
    HcalCoderDb coder(*channelCoder, *shape);
    CaloSamples samples;
    coder.adc2fC(digi, samples);

    for (int ts = 0; ts < samples.size(); ++ts) {
      it->second->Fill(digi[ts].capid(), samples[ts]);
    }
  }
}

void HcalPedestalCalibDQMAnalyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  HcalDbService const& dbService = iSetup.getData(tokenDbService_);

  edm::Handle<QIE11DigiCollection> qie11Digis;
  iEvent.getByToken(tokenQIE11_, qie11Digis);
  edm::Handle<QIE10DigiCollection> qie10Digis;
  iEvent.getByToken(tokenQIE10_, qie10Digis);
  edm::Handle<HODigiCollection> hoDigis;
  iEvent.getByToken(tokenHO_, hoDigis);

  fillHBHE(*qie11Digis, HcalBarrel, dbService);
  fillHBHE(*qie11Digis, HcalEndcap, dbService);
  fillHF(*qie10Digis, dbService);
  fillHO(*hoDigis, dbService);
}

DEFINE_FWK_MODULE(HcalPedestalCalibDQMAnalyzer);
