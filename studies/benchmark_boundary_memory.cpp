#include <sys/resource.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "OversampledHistogram.h"

namespace {
  long PeakRssKiB() {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
      throw std::runtime_error("getrusage failed");
    }
    return usage.ru_maxrss;
  }
}  // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 7) {
      throw std::invalid_argument("Usage: benchmark_boundary_memory EVENTS FACTOR BINS SLOTS RANGE_ROWS OCCUPANCY");
    }
    const unsigned long events = std::stoul(argv[1]);
    const int factor = std::stoi(argv[2]);
    const int bins = std::stoi(argv[3]);
    const unsigned int slots = std::stoul(argv[4]);
    const unsigned long rangeRows = std::stoul(argv[5]);
    const int occupancy = std::stoi(argv[6]);
    if (!events || factor < 1 || bins < 1 || !slots || !rangeRows || occupancy < 1 || occupancy > bins) {
      throw std::invalid_argument("Invalid benchmark dimensions");
    }
    ROOT::EnableImplicitMT(slots);
    {
      RangeAwareOversampledHistogram1D warm(1, "warm", "", 1, 0, 1);
      warm.Exec(0, 0UL, 0.5, 0UL, 1UL);
      warm.Finalize();
    }
    ROOT::VecOps::RVec<double> values;
    for (int bin = 0; bin < occupancy; ++bin) {
      values.push_back(bin + 0.5);
    }
    const auto rssBefore = PeakRssKiB();
    RangeAwareOversampledHistogram1D action(factor, "memory", "", bins, 0, bins);
    const auto rssSlots = PeakRssKiB();
    const auto start = std::chrono::steady_clock::now();
    const auto rows = events * factor;
    // Deterministic round-robin range assignment; serialized calls isolate storage costs.
    for (unsigned long row = 0; row < rows; ++row) {
      const auto range = row / rangeRows;
      action.Exec(range % slots, row / factor, values, range * rangeRows, std::min((range + 1) * rangeRows, rows));
    }
    const auto filled = std::chrono::steady_clock::now();
    const auto rssFilled = PeakRssKiB();
    action.Finalize();
    const auto finalized = std::chrono::steady_clock::now();
    for (int bin = 1; bin <= occupancy; ++bin) {
      if (action.GetResultPtr()->GetBinContent(bin) != events ||
          std::abs(action.GetResultPtr()->GetBinError(bin) - std::sqrt(events)) > 1e-8) {
        throw std::runtime_error("Benchmark changed the expected event-level contents or errors");
      }
    }
    if (action.GetResultPtr()->GetEntries() != events * occupancy) {
      throw std::runtime_error("Benchmark changed the expected fill count");
    }
    std::cout << "RESULT," << events << ',' << factor << ',' << bins << ',' << slots << ',' << rangeRows << ','
              << occupancy << ',' << std::chrono::duration<double>(filled - start).count() << ','
              << std::chrono::duration<double>(finalized - filled).count() << ',' << rssBefore << ',' << rssSlots << ','
              << rssFilled << ',' << PeakRssKiB() << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
