// system includes
#include <algorithm>
#include <cmath>
#include <format>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// ROOT includes
#include "TF1.h"
#include "TH1D.h"
#include "TH2.h"

// user includes
#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

// ---------------------------------------------------------------------------
// Robust E/p vs eta from the E/p vs eta_SC 2D histograms of ScoutingDileptonMonitor. In each eta_SC bin:
//  - median and 16%/84% quantiles;
//  - shorth: midpoint of the shortest interval containing 50% of the entries, and sigma_eff: half-width
//    of the shortest interval containing 68.3%. They follow the peak, while the median is pulled towards
//    the one-sided bremsstrahlung tail in the endcaps;
//  - optionally, mean and width of an iterative gaussian fit of the core, seeded by shorth and sigma_eff.
// ---------------------------------------------------------------------------
namespace {

  // Cumulative distribution of E/p in one eta slice, piecewise linear between the y bin edges (uniform
  // density inside each bin). Underflow and overflow are point masses at the axis edges, so the quantiles
  // are those of all the entries; a quantile falling in them is clamped to the axis range.
  class SliceCDF {
  public:
    SliceCDF(const TH2& h, int ix) {
      const TAxis* yAxis = h.GetYaxis();
      const int ny = h.GetNbinsY();
      edges_.reserve(ny + 1);
      cumulative_.reserve(ny + 1);
      double sum = h.GetBinContent(ix, 0);
      edges_.push_back(yAxis->GetBinLowEdge(1));
      cumulative_.push_back(sum);
      for (int iy = 1; iy <= ny; ++iy) {
        sum += h.GetBinContent(ix, iy);
        edges_.push_back(yAxis->GetBinUpEdge(iy));
        cumulative_.push_back(sum);
      }
      total_ = sum + h.GetBinContent(ix, ny + 1);
    }

    double total() const { return total_; }

    // y below which a fraction p of the entries lies. Across empty bins the quantile function jumps:
    // the default returns the lower end of the jump, upperLimit = true the upper one.
    double quantile(double p, bool upperLimit = false) const {
      const double target = p * total_;
      const auto it = upperLimit ? std::upper_bound(cumulative_.begin(), cumulative_.end(), target)
                                 : std::lower_bound(cumulative_.begin(), cumulative_.end(), target);
      if (it == cumulative_.begin())
        return edges_.front();
      if (it == cumulative_.end())
        return edges_.back();
      const size_t k = std::distance(cumulative_.begin(), it);
      const double c0 = cumulative_[k - 1], c1 = cumulative_[k];
      return edges_[k - 1] + (edges_[k] - edges_[k - 1]) * (target - c0) / (c1 - c0);
    }

    // Shortest interval containing a fraction f of the entries. Its width is piecewise linear in the
    // fraction p below it, so the minimum is where one of the two ends sits on a bin edge.
    std::pair<double, double> shortestInterval(double f) const {
      std::pair<double, double> best{edges_.front(), edges_.back()};
      auto test = [&](double p) {
        p = std::clamp(p, 0., 1. - f);
        const double lo = quantile(p, true);
        const double hi = quantile(p + f);
        if (hi - lo < best.second - best.first)
          best = {lo, hi};
      };
      for (const double c : cumulative_) {
        test(c / total_);
        test(c / total_ - f);
      }
      return best;
    }

  private:
    std::vector<double> edges_;
    std::vector<double> cumulative_;  // entries below each edge, underflow included
    double total_{0.};
  };

}  // namespace

class ScoutingEOverPHarvester : public DQMEDHarvester {
public:
  explicit ScoutingEOverPHarvester(const edm::ParameterSet&);
  ~ScoutingEOverPHarvester() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override;

private:
  bool coreFit(TH1& proj, TF1& gaus, double& mean, double& meanError, double& sigma) const;

  const std::string folder_;
  const std::vector<std::string> inputHistos_;
  const double minEntries_;
  const bool doCoreFit_;
  const double minEntriesCoreFit_;
  const double coreFitNSigma_;
  const unsigned int coreFitIterations_;
};

ScoutingEOverPHarvester::ScoutingEOverPHarvester(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      inputHistos_(iConfig.getParameter<std::vector<std::string>>("inputHistos")),
      minEntries_(iConfig.getParameter<double>("minEntries")),
      doCoreFit_(iConfig.getParameter<bool>("doCoreFit")),
      minEntriesCoreFit_(iConfig.getParameter<double>("minEntriesCoreFit")),
      coreFitNSigma_(iConfig.getParameter<double>("coreFitNSigma")),
      coreFitIterations_(iConfig.getParameter<unsigned int>("coreFitIterations")) {}

// Iterative gaussian fit within mean +- coreFitNSigma sigma; mean and sigma are the seeds on input
bool ScoutingEOverPHarvester::coreFit(TH1& proj, TF1& gaus, double& mean, double& meanError, double& sigma) const {
  for (unsigned int i = 0; i < coreFitIterations_; ++i) {
    if (!(sigma > 0.))
      return false;
    gaus.SetRange(mean - coreFitNSigma_ * sigma, mean + coreFitNSigma_ * sigma);
    gaus.SetParameters(proj.GetBinContent(proj.FindBin(mean)), mean, sigma);
    if (static_cast<int>(proj.Fit(&gaus, "QNR")) != 0)
      return false;
    mean = gaus.GetParameter(1);
    meanError = gaus.GetParError(1);
    sigma = std::abs(gaus.GetParameter(2));
  }
  return mean > proj.GetXaxis()->GetXmin() && mean < proj.GetXaxis()->GetXmax();
}

void ScoutingEOverPHarvester::dqmEndJob(DQMStore::IBooker& ibooker, DQMStore::IGetter& igetter) {
  ibooker.setCurrentFolder(folder_);
  TF1 gaus("scoutingEOverPCoreGaus", "gaus", 0., 1.);

  for (const auto& inputName : inputHistos_) {
    const MonitorElement* me = igetter.get(folder_ + "/" + inputName);
    if (me == nullptr || me->getTH2F() == nullptr) {
      edm::LogInfo("ScoutingEOverPHarvester") << "Missing 2D histogram " << folder_ << "/" << inputName << ", skipping";
      continue;
    }
    const TH2F& h2 = *me->getTH2F();

    // "<name>_EoverP_vsEta2D" -> "<name>_EoverP_vsEta_<estimator>"
    std::string base = inputName;
    if (base.size() > 2 && base.ends_with("2D"))
      base.erase(base.size() - 2);

    const TAxis* xAxis = h2.GetXaxis();
    const int nx = xAxis->GetNbins();
    const std::string title = h2.GetTitle();
    auto book = [&](const std::string& suffix, const std::string& what) {
      return ibooker.book1D(base + suffix,
                            title + ", " + what + ";" + xAxis->GetTitle() + ";" + what + " of E_{SC}/p_{trk}",
                            nx,
                            xAxis->GetXmin(),
                            xAxis->GetXmax());
    };
    MonitorElement* median = book("_median", "median");
    MonitorElement* q16 = book("_q16", "16% quantile");
    MonitorElement* q84 = book("_q84", "84% quantile");
    MonitorElement* shorth = book("_shorth", "shortest-50% midpoint");
    MonitorElement* sigmaEff = book("_sigmaEff", "#sigma_{eff}");
    MonitorElement* coreMean = doCoreFit_ ? book("_coreMean", "core gaussian mean") : nullptr;
    MonitorElement* coreSigma = doCoreFit_ ? book("_coreSigma", "core gaussian #sigma") : nullptr;

    for (int ix = 1; ix <= nx; ++ix) {
      const SliceCDF cdf(h2, ix);
      const double total = cdf.total();
      if (total < minEntries_ || total <= 0.)
        continue;

      const double lo = cdf.quantile(0.16);
      const double med = cdf.quantile(0.50);
      const double hi = cdf.quantile(0.84);
      // uncertainty of the median for a gaussian core: sqrt(pi/2) sigma / sqrt(N), sigma from the quantiles
      median->setBinContent(ix, med);
      median->setBinError(ix, std::sqrt(M_PI / 2.) * 0.5 * (hi - lo) / std::sqrt(total));
      q16->setBinContent(ix, lo);
      q84->setBinContent(ix, hi);

      const auto [lo50, hi50] = cdf.shortestInterval(0.5);
      const auto [lo68, hi68] = cdf.shortestInterval(0.683);
      const double mid = 0.5 * (lo50 + hi50);
      const double sEff = 0.5 * (hi68 - lo68);
      shorth->setBinContent(ix, mid);
      sigmaEff->setBinContent(ix, sEff);

      if (!doCoreFit_ || total < minEntriesCoreFit_)
        continue;
      std::unique_ptr<TH1D> proj(h2.ProjectionY(std::format("{}_slice{}", base, ix).c_str(), ix, ix));
      proj->SetDirectory(nullptr);
      double mean = mid, meanError = 0., sigma = sEff;
      if (coreFit(*proj, gaus, mean, meanError, sigma)) {
        coreMean->setBinContent(ix, mean);
        coreMean->setBinError(ix, meanError);
        coreSigma->setBinContent(ix, sigma);
      } else {
        edm::LogInfo("ScoutingEOverPHarvester") << "Core fit failed for " << inputName << " eta bin " << ix;
      }
    }
  }
}

void ScoutingEOverPHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("folder", "HLT/ScoutingOffline/DiLepton");
  desc.add<std::vector<std::string>>("inputHistos", {"electronsZ_EoverP_vsEta2D", "electronsJpsi_EoverP_vsEta2D"})
      ->setComment("E/p vs eta_SC 2D histograms booked by ScoutingDileptonMonitor in folder");
  desc.add<double>("minEntries", 20.)->setComment("eta bins with fewer entries are left empty");
  desc.add<bool>("doCoreFit", true)->setComment("iterative gaussian fit of the E/p core, seeded by shorth and sigma_eff");
  desc.add<double>("minEntriesCoreFit", 200.)->setComment("eta bins with fewer entries are not fitted");
  desc.add<double>("coreFitNSigma", 1.5)->setComment("fit range: mean +- coreFitNSigma sigma");
  desc.add<unsigned int>("coreFitIterations", 3);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(ScoutingEOverPHarvester);
