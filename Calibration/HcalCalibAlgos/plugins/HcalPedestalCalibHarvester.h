#ifndef Calibration_HcalCalibAlgos_HcalPedestalCalibHarvester_h
#define Calibration_HcalCalibAlgos_HcalPedestalCalibHarvester_h

#include "DQMServices/Core/interface/DQMEDHarvester.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

class HcalTopology;
class HcalRecNumberingRecord;

class HcalPedestalCalibHarvester : public DQMEDHarvester {
public:
  explicit HcalPedestalCalibHarvester(const edm::ParameterSet&);
  ~HcalPedestalCalibHarvester() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

protected:
  void dqmEndRun(DQMStore::IBooker&, DQMStore::IGetter&, edm::Run const&, edm::EventSetup const&) override;
  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override {}

private:
  const double minEntriesPerCapId_;
  const std::string dqmDir_;
  const std::string recordLabel_;
  edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> tokenTopology_;
};

#endif
