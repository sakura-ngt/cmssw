#ifndef Calibration_HcalCalibAlgos_HcalPedestalCalibDQMAnalyzer_h
#define Calibration_HcalCalibAlgos_HcalPedestalCalibDQMAnalyzer_h

#include <unordered_map>

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

class MonitorElement;

class HcalPedestalCalibDQMAnalyzer : public DQMEDAnalyzer {
public:
  explicit HcalPedestalCalibDQMAnalyzer(const edm::ParameterSet&);
  ~HcalPedestalCalibDQMAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

protected:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(edm::Event const&, edm::EventSetup const&) override;

private:
  // Books one capId-vs-fC TH2 per channel for the given subdet, using the
  // run-level DetId list product (mirrors HcalDigiSortedTableProducer::beginRun).
  void bookSubdet(DQMStore::IBooker&,
                  edm::Run const&,
                  edm::EDGetTokenT<std::vector<HcalDetId>> const&,
                  HcalSubdetector,
                  std::string const& folderSuffix);

  void fillHBHE(QIE11DigiCollection const&, HcalSubdetector, HcalDbService const&);
  void fillHF(QIE10DigiCollection const&, HcalDbService const&);
  void fillHO(HODigiCollection const&, HcalDbService const&);

  std::unordered_map<uint32_t, MonitorElement*> pedestalByCapId_;

  edm::EDGetTokenT<std::vector<HcalDetId>> tokenHBDetIdList_;
  edm::EDGetTokenT<std::vector<HcalDetId>> tokenHEDetIdList_;
  edm::EDGetTokenT<std::vector<HcalDetId>> tokenHFDetIdList_;
  edm::EDGetTokenT<std::vector<HcalDetId>> tokenHODetIdList_;

  edm::EDGetTokenT<QIE11DigiCollection> tokenQIE11_;
  edm::EDGetTokenT<QIE10DigiCollection> tokenQIE10_;
  edm::EDGetTokenT<HODigiCollection> tokenHO_;

  edm::ESGetToken<HcalDbService, HcalDbRecord> tokenDbService_;

  const int fcAxisBins_;
  const double fcAxisMax_;
};

#endif
