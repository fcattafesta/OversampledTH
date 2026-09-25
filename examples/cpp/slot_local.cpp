#include <ROOT/RDataFrame.hxx>
#include <iostream>

#include "OversampledHistogram.h"

int main() {
  ROOT::EnableImplicitMT(2);
  ROOT::RDataFrame frame(100);
  auto rdf = frame.Define("genEvent", "static_cast<unsigned long>(rdfentry_ / 2)").Define("value", "0.5");

  SlotLocalOversampledHistogram1D helper(2, "slot_local", "Example", 2, 0, 2);
  auto histogram = rdf.Book<unsigned long, double>(helper, {"genEvent", "value"});
  const auto content = histogram->GetBinContent(1);  // Triggers the event loop and diagnostic.
  std::cout << "Content: " << content << ", error: " << histogram->GetBinError(1) << '\n';
  // Check the split count printed at finalization before relying on its errors.
}
