#include <ROOT/RDataFrame.hxx>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "OversampledHistogram.h"

namespace {

  void Check(bool condition, const std::string &message) {
    if (!condition) {
      throw std::runtime_error(message);
    }
  }

  void Near(double actual, double expected, const std::string &message) {
    Check(std::abs(actual - expected) < 1e-10,
          message + ": got " + std::to_string(actual) + ", expected " + std::to_string(expected));
  }

  void TestSequential() {
    Check(!ROOT::IsImplicitMTEnabled(), "sequential test requires implicit MT to be disabled");
    SequentialOversampledHistogram1D action(2, "sequential", "", 2, 0, 2);
    action.Exec(0, 1UL, 0.5, 2.0);
    action.Exec(0, 1UL, 0.5, 4.0);
    action.Exec(0, 2UL, 0.5, -1.0);
    action.Exec(0, 2UL, ROOT::VecOps::RVec<double>{1.5, 1.5});
    action.Finalize();

    const auto result = action.GetResultPtr();
    Near(result->GetBinContent(1), 2.5, "sequential content");
    Near(result->GetBinError(1), std::sqrt(9.25), "sequential correlated error");
    Near(result->GetBinContent(2), 1.0, "sequential vector content");
    Near(result->GetBinError(2), 1.0, "sequential vector error");
    action.Finalize();
    Near(result->GetBinContent(1), 2.5, "sequential finalize is idempotent");

    SequentialOversampledHistogram1D edge(2, "edge", "", 2, 0, 2);
    edge.Exec(0, 1UL, -1.0, 2.0);
    edge.Exec(0, 1UL, 3.0, 4.0);
    edge.Finalize();
    Near(edge.GetResultPtr()->GetBinContent(0), 1.0, "underflow content");
    Near(edge.GetResultPtr()->GetBinError(0), 1.0, "underflow error");
    Near(edge.GetResultPtr()->GetBinContent(3), 2.0, "overflow content");
    Near(edge.GetResultPtr()->GetBinError(3), 2.0, "overflow error");

    SequentialOversampledHistogram1D copied(action);
    Near(copied.GetResultPtr()->GetBinContent(1), 0.0, "copy starts empty");
    copied = action;
    Near(copied.GetResultPtr()->GetBinContent(1), 0.0, "copy assignment starts empty");
    bool rejected = false;
    try {
      SequentialOversampledHistogram1D invalid(0, "invalid", "", 2, 0, 2);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    Check(rejected, "zero oversampling factor must be rejected");
  }

  void TestRDataFrameBooking() {
    ROOT::RDataFrame df(4);
    auto source = df.Define("event", "rdfentry_ < 2 ? 1UL : 2UL").Define("value", "0.5");
    SequentialOversampledHistogram1D action(2, "booked", "", 2, 0, 2);
    auto result = source.Book<unsigned long, double>(action, {"event", "value"});
    Near(result->GetBinContent(1), 2.0, "RDataFrame booking content");
    Near(result->GetBinError(1), std::sqrt(2.0), "RDataFrame booking error");
  }

  void TestMultiThreadedActions() {
    ROOT::EnableImplicitMT(2);
    bool rejected = false;
    try {
      SequentialOversampledHistogram1D invalid(2, "invalid_mt", "", 2, 0, 2);
    } catch (const std::runtime_error &) {
      rejected = true;
    }
    Check(rejected, "sequential action must reject implicit MT");

    SlotLocalOversampledHistogram1D local(2, "local", "", 2, 0, 2);
    local.Exec(0, 1UL, 0.5, 2.0);
    local.Exec(1, 1UL, 0.5, 4.0);
    local.Finalize();
    Near(local.GetResultPtr()->GetBinContent(1), 3.0, "slot-local split content");
    Near(local.GetResultPtr()->GetBinError(1), std::sqrt(5.0), "slot-local split error");
    SlotLocalOversampledHistogram1D localCopy(local);
    localCopy = local;
    Near(localCopy.GetResultPtr()->GetBinContent(1), 0.0, "slot-local copy assignment");

    SlotLocalOversampledHistogram1D revisited(2, "revisited", "", 2, 0, 2);
    revisited.Exec(0, 1UL, 0.5, 2.0);
    revisited.Exec(0, 2UL, 0.5, 1.0);
    revisited.Exec(0, 1UL, 0.5, 4.0);
    std::ostringstream diagnostics;
    auto *oldBuffer = std::cout.rdbuf(diagnostics.rdbuf());
    revisited.Finalize();
    std::cout.rdbuf(oldBuffer);
    Check(diagnostics.str().find("split across slots=0, revisited within slots=1") != std::string::npos,
          "same-slot revisit must be reported");
    Near(revisited.GetResultPtr()->GetBinError(1), std::sqrt(5.25), "same-slot revisit error");

    RangeAwareOversampledHistogram1D range(2, "range", "", 2, 0, 2);
    range.Exec(0, 1UL, 0.5, 0UL, 2UL, 2.0);
    range.Exec(1, 1UL, 0.5, 2UL, 6UL, 4.0);
    range.Exec(1, 2UL, 0.5, 2UL, 6UL, 2.0);
    range.Exec(1, 3UL, 0.5, 2UL, 6UL, 1.0);
    range.Finalize();
    Near(range.GetResultPtr()->GetBinContent(1), 4.5, "range-aware split content");
    Near(range.GetResultPtr()->GetBinError(1), std::sqrt(10.25), "range-aware split error");
    RangeAwareOversampledHistogram1D reusedSlot(2, "reused_slot", "", 2, 0, 2);
    reusedSlot.Exec(0, 1UL, 0.5, 0UL, 2UL, 2.0);
    reusedSlot.Exec(0, 1UL, 0.5, 2UL, 4UL, 4.0);
    reusedSlot.Finalize();
    Near(reusedSlot.GetResultPtr()->GetBinError(1), 3.0, "range split on one reused slot");
    RangeAwareOversampledHistogram1D rangeCopy(range);
    rangeCopy = range;
    Near(rangeCopy.GetResultPtr()->GetBinContent(1), 0.0, "range-aware copy assignment");
    ROOT::DisableImplicitMT();
  }

}  // namespace

int main() {
  try {
    TestSequential();
    TestRDataFrameBooking();
    TestMultiThreadedActions();
    std::cout << "All histogram tests passed\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
