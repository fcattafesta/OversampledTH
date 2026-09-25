#include <sys/resource.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "OversampledHistogram.h"
namespace {

  constexpr int kFactor = 3;
  constexpr int kRangeRows = 257;

  long PeakRssKiB() {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
      throw std::runtime_error("getrusage failed");
    }
    return usage.ru_maxrss;
  }

  template <typename Action, typename Fill>
  double Run(Action &action, int events, int threads, Fill fill) {
    const auto start = std::chrono::steady_clock::now();
    if (threads == 1) {
      fill(0, 0, events);
    } else {
      std::vector<std::thread> workers;
      workers.reserve(threads);
      for (int slot = 0; slot < threads; ++slot) {
        const int begin = events * slot / threads;
        const int end = events * (slot + 1) / threads;
        workers.emplace_back(fill, slot, begin, end);
      }
      for (auto &worker : workers) {
        worker.join();
      }
    }
    action.Finalize();
    const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (std::abs(action.GetResultPtr()->GetBinContent(1) - events) > 1e-8) {
      throw std::runtime_error("Benchmark result failed its event-count check");
    }
    return seconds;
  }

}  // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 4) {
      throw std::invalid_argument("Usage: benchmark_actions MODE EVENTS THREADS");
    }
    const std::string mode = argv[1];
    const int events = std::stoi(argv[2]);
    const int threads = std::stoi(argv[3]);
    if (events < 1 || threads < 1 || (mode == "sequential" && threads != 1)) {
      throw std::invalid_argument("Invalid event or thread count");
    }
    if (mode != "sequential" && mode != "slot-local" && mode != "range-aware") {
      throw std::invalid_argument("Mode must be sequential, slot-local, or range-aware");
    }
    if (threads > 1) {
      ROOT::EnableImplicitMT(threads);
    }

    // Initialize ROOT's histogram machinery before measuring RSS growth.
    if (mode == "sequential") {
      SequentialOversampledHistogram1D warm(1, "warm_seq", "", 1, 0, 1);
      warm.Exec(0, 0UL, 0.5);
      warm.Finalize();
    } else if (mode == "slot-local") {
      SlotLocalOversampledHistogram1D warm(1, "warm_slot", "", 1, 0, 1);
      warm.Exec(0, 0UL, 0.5);
      warm.Finalize();
    } else {
      RangeAwareOversampledHistogram1D warm(1, "warm_range", "", 1, 0, 1);
      warm.Exec(0, 0UL, 0.5, 0UL, 1UL);
      warm.Finalize();
    }

    const auto rssBefore = PeakRssKiB();
    double seconds = 0.0;
    if (mode == "sequential") {
      SequentialOversampledHistogram1D action(kFactor, "bench_seq", "", 1, 0, 1);
      seconds = Run(action, events, 1, [&](int, int begin, int end) {
        for (int event = begin; event < end; ++event) {
          for (int replica = 0; replica < kFactor; ++replica) {
            action.Exec(0, static_cast<unsigned long>(event), 0.5);
          }
        }
      });
    } else if (mode == "slot-local") {
      SlotLocalOversampledHistogram1D action(kFactor, "bench_slot", "", 1, 0, 1);
      seconds = Run(action, events, threads, [&](int slot, int begin, int end) {
        for (int event = begin; event < end; ++event) {
          for (int replica = 0; replica < kFactor; ++replica) {
            action.Exec(slot, static_cast<unsigned long>(event), 0.5);
          }
        }
      });
    } else if (mode == "range-aware") {
      RangeAwareOversampledHistogram1D action(kFactor, "bench_range", "", 1, 0, 1);
      seconds = Run(action, events, threads, [&](int slot, int begin, int end) {
        for (int event = begin; event < end; ++event) {
          for (int replica = 0; replica < kFactor; ++replica) {
            const int localRow = (event - begin) * kFactor + replica;
            const auto rangeBegin = static_cast<unsigned long>(begin * kFactor + (localRow / kRangeRows) * kRangeRows);
            const auto rangeEnd = std::min<unsigned long>(rangeBegin + kRangeRows, end * kFactor);
            action.Exec(slot, static_cast<unsigned long>(event), 0.5, rangeBegin, rangeEnd);
          }
        }
      });
    }

    const long rssGrowth = PeakRssKiB() - rssBefore;
    std::cout << "RESULT," << mode << ',' << events << ',' << threads << ',' << seconds << ','
              << (events * kFactor / seconds) << ',' << rssGrowth << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
