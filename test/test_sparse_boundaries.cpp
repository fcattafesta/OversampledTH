#include <sys/resource.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>

#include "OversampledHistogram.h"

namespace {

  template <typename TH>
  void TestDenseEquivalence() {
    constexpr int bins = 64;
    constexpr int slots = 8;
    constexpr int factor = 12;
    RangeAwareOversampledHistogram<TH> action(factor, "sparse", "", bins, 0, bins);
    TH model("dense_model", "", bins, 0, bins);
    model.SetDirectory(nullptr);
    model.Sumw2();
    using Map = std::unordered_map<unsigned long, std::unique_ptr<TH>>;
    std::vector<Map> boundaries(slots);
    unsigned long row = 0;
    for (unsigned long event = 0; event < 4; ++event) {
      for (int piece = 0; piece < factor; ++piece, ++row) {
        ROOT::VecOps::RVec<double> values{-1.0, 0.5, 1.5, 70.0};
        double weight = (piece % 3 + 1) * 0.1;
        if (event == 1) {
          weight = piece % 2 == 0 ? 1.0 : -1.0;
        } else if (event == 2 && piece % 2 == 0) {
          values.clear();
        } else if (event == 3) {
          values.clear();
          for (int bin = 0; bin < bins; ++bin) {
            values.push_back(bin + 0.5);
          }
        }
        const auto slot = static_cast<unsigned int>(row % slots);
        action.Exec(slot, event, values, row, row + 1, weight);
        auto partial = oversampled_detail::CloneEmpty(model);
        for (const auto value : values) {
          partial->Fill(value, weight);
        }
        auto &target = boundaries[slot][event];
        if (!target) {
          target = oversampled_detail::CloneEmpty(model);
        }
        target->Add(partial.get());
      }
    }
    Map merged;
    for (auto &slot : boundaries) {
      for (auto &[id, hist] : slot) {
        auto &target = merged[id];
        if (!target) {
          target = oversampled_detail::CloneEmpty(model);
        }
        target->Add(hist.get());
      }
    }
    auto reference = oversampled_detail::CloneEmpty(model);
    for (const auto &[id, hist] : merged) {
      (void)id;
      oversampled_detail::FillEvent(*reference, *hist, factor);
    }
    action.Finalize();
    const auto check = [&] {
      const double tolerance = std::is_same_v<TH, TH1F> ? 1e-5 : 1e-11;
      for (int bin = 0; bin <= bins + 1; ++bin) {
        if (std::abs(action.GetResultPtr()->GetBinContent(bin) - reference->GetBinContent(bin)) > tolerance ||
            std::abs(action.GetResultPtr()->GetBinError(bin) - reference->GetBinError(bin)) > tolerance) {
          throw std::runtime_error("Sparse boundary result differs from dense reference");
        }
      }
      if (action.GetResultPtr()->GetEntries() != reference->GetEntries()) {
        throw std::runtime_error("Sparse boundary result changed the event-level fill count");
      }
    };
    check();
    action.Finalize();
    check();
  }

  long PeakRssKiB() {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
      throw std::runtime_error("getrusage failed");
    }
    return usage.ru_maxrss;
  }

  void TestLargeBinMemory() {
    constexpr int bins = 10000;
    constexpr int events = 1000;
    constexpr int factor = 9;
    constexpr int rangeRows = 7;
    const auto before = PeakRssKiB();
    RangeAwareOversampledHistogram1D action(factor, "large_bins", "", bins, 0, bins);
    for (unsigned long row = 0; row < events * factor; ++row) {
      const auto range = row / rangeRows;
      action.Exec(range % 8, row / factor, 0.5, range * rangeRows, (range + 1) * rangeRows);
    }
    action.Finalize();
    const auto growth = PeakRssKiB() - before;
    if (action.GetResultPtr()->GetBinContent(1) != events ||
        std::abs(action.GetResultPtr()->GetBinError(1) - std::sqrt(events)) > 1e-10) {
      throw std::runtime_error("Large-bin stress result changed contents or Sumw2 errors");
    }
    if (growth > 64 * 1024) {
      throw std::runtime_error("Large-bin boundary storage exceeded the 64 MiB RSS growth budget");
    }
    std::cout << "Large-bin boundary peak RSS growth: " << growth << " KiB\n";
  }

}  // namespace

int main() {
  try {
    ROOT::EnableImplicitMT(8);
    TestDenseEquivalence<TH1D>();
    TestDenseEquivalence<TH1F>();
    TestLargeBinMemory();
    ROOT::DisableImplicitMT();
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
