#include <sys/resource.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "OversampledHistogram.h"

namespace {

  constexpr int kEvents = 50000;
  constexpr int kFactor = 3;
  constexpr int kRangeRows = 257;
  constexpr long kMemoryBudgetKiB = 128 * 1024;

  long PeakRssKiB() {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
      throw std::runtime_error("getrusage failed");
    }
    return usage.ru_maxrss;  // KiB on Linux.
  }

  void CheckHistogram(const TH1D &hist, const char *label) {
    if (std::abs(hist.GetBinContent(1) - kEvents) > 1e-8 || std::abs(hist.GetBinError(1) - std::sqrt(kEvents)) > 1e-8 ||
        std::abs(hist.GetEntries() - kEvents) > 1e-8) {
      throw std::runtime_error(std::string(label) + " lost or duplicated generated events");
    }
  }

}  // namespace

int main() {
  try {
    ROOT::EnableImplicitMT(2);
    const auto before = PeakRssKiB();
    SlotLocalOversampledHistogram1D slotLocal(kFactor, "slot_stress", "", 1, 0, 1);
    RangeAwareOversampledHistogram1D rangeAware(kFactor, "range_stress", "", 1, 0, 1);

    int row = 0;
    for (int event = 0; event < kEvents; ++event) {
      for (int replica = 0; replica < kFactor; ++replica, ++row) {
        slotLocal.Exec(event % 2, static_cast<unsigned long>(event), 0.5);
        const auto rangeBegin = static_cast<unsigned long>((row / kRangeRows) * kRangeRows);
        const auto rangeEnd = std::min<unsigned long>(rangeBegin + kRangeRows, kEvents * kFactor);
        rangeAware.Exec((row / kRangeRows) % 2, static_cast<unsigned long>(event), 0.5, rangeBegin, rangeEnd);
      }
    }
    slotLocal.Finalize();
    rangeAware.Finalize();
    CheckHistogram(*slotLocal.GetResultPtr(), "slot-local");
    CheckHistogram(*rangeAware.GetResultPtr(), "range-aware");

    const auto growth = PeakRssKiB() - before;
    std::cout << "Processed " << kEvents << " generated events; peak RSS growth " << growth << " KiB\n";
    if (growth > kMemoryBudgetKiB) {
      throw std::runtime_error("Peak RSS growth exceeded the 128 MiB stress budget");
    }
    ROOT::DisableImplicitMT();
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
