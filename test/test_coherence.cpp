#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <thread>

#include "OversampledHistogram.h"

namespace {

  constexpr int kEvents = 128;
  constexpr int kFactor = 3;
  constexpr int kBins = 4;

  void Compare(const TH1D &actual, const TH1D &expected, const char *label) {
    for (int bin = 0; bin <= kBins + 1; ++bin) {
      if (std::abs(actual.GetBinContent(bin) - expected.GetBinContent(bin)) > 1e-9 ||
          std::abs(actual.GetBinError(bin) - expected.GetBinError(bin)) > 1e-9) {
        throw std::runtime_error(std::string(label) + " differs in bin " + std::to_string(bin));
      }
    }
  }

  double Weight(int event, int replica) { return (event % 7 == 0 ? -0.5 : 1.0) + replica * 0.25; }

  double Value(int event, int replica) { return static_cast<double>((event + replica) % kBins) + 0.5; }

}  // namespace

int main() {
  try {
    TH1D expected("expected", "", kBins, 0, kBins);
    expected.SetDirectory(nullptr);
    expected.Sumw2();

    SequentialOversampledHistogram1D sequential(kFactor, "sequential", "", kBins, 0, kBins);
    for (int event = 0; event < kEvents; ++event) {
      std::array<double, kBins> totals{};
      for (int replica = 0; replica < kFactor; ++replica) {
        const auto value = Value(event, replica);
        const auto weight = Weight(event, replica);
        totals[static_cast<int>(value)] += weight;
        sequential.Exec(0, static_cast<unsigned long>(event), value, weight);
      }
      for (int bin = 0; bin < kBins; ++bin) {
        if (totals[bin] != 0.0) {
          expected.Fill(bin + 0.5, totals[bin] / kFactor);
        }
      }
    }
    sequential.Finalize();
    Compare(*sequential.GetResultPtr(), expected, "sequential");

    ROOT::EnableImplicitMT(2);
    SlotLocalOversampledHistogram1D slotLocal(kFactor, "slot", "", kBins, 0, kBins);
    RangeAwareOversampledHistogram1D rangeAware(kFactor, "range", "", kBins, 0, kBins);
    int row = 0;
    for (int event = 0; event < kEvents; ++event) {
      for (int replica = 0; replica < kFactor; ++replica, ++row) {
        const auto value = Value(event, replica);
        const auto weight = Weight(event, replica);
        slotLocal.Exec(event % 2, static_cast<unsigned long>(event), value, weight);
        const auto rangeBegin = static_cast<unsigned long>((row / 7) * 7);
        const auto rangeEnd = std::min<unsigned long>(rangeBegin + 7, kEvents * kFactor);
        rangeAware.Exec((row / 7) % 2, static_cast<unsigned long>(event), value, rangeBegin, rangeEnd, weight);
      }
    }
    slotLocal.Finalize();
    rangeAware.Finalize();
    Compare(*slotLocal.GetResultPtr(), expected, "slot-local");
    Compare(*rangeAware.GetResultPtr(), expected, "range-aware");

    SlotLocalOversampledHistogram1D parallelSlot(kFactor, "parallel_slot", "", 1, 0, 1);
    std::thread local0([&] {
      for (int event = 0; event < kEvents / 2; ++event) {
        for (int replica = 0; replica < kFactor; ++replica) {
          parallelSlot.Exec(0, static_cast<unsigned long>(event), 0.5);
        }
      }
    });
    std::thread local1([&] {
      for (int event = kEvents / 2; event < kEvents; ++event) {
        for (int replica = 0; replica < kFactor; ++replica) {
          parallelSlot.Exec(1, static_cast<unsigned long>(event), 0.5);
        }
      }
    });
    local0.join();
    local1.join();
    parallelSlot.Finalize();

    RangeAwareOversampledHistogram1D parallelRange(kFactor, "parallel_range", "", 1, 0, 1);
    constexpr int splitRow = (kEvents / 2) * kFactor + 1;
    std::thread range0([&] {
      for (int index = 0; index < splitRow; ++index) {
        parallelRange.Exec(0, static_cast<unsigned long>(index / kFactor), 0.5, 0UL, splitRow);
      }
    });
    std::thread range1([&] {
      for (int index = splitRow; index < kEvents * kFactor; ++index) {
        parallelRange.Exec(1, static_cast<unsigned long>(index / kFactor), 0.5, splitRow, kEvents * kFactor);
      }
    });
    range0.join();
    range1.join();
    parallelRange.Finalize();
    for (const auto *hist : {parallelSlot.GetResultPtr().get(), parallelRange.GetResultPtr().get()}) {
      if (std::abs(hist->GetBinContent(1) - kEvents) > 1e-9 ||
          std::abs(hist->GetBinError(1) - std::sqrt(kEvents)) > 1e-9 || std::abs(hist->GetEntries() - kEvents) > 1e-9) {
        throw std::runtime_error("Parallel processing lost or duplicated an event");
      }
    }
    ROOT::DisableImplicitMT();
    std::cout << "Coherence checks passed for " << kEvents << " generated events\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
