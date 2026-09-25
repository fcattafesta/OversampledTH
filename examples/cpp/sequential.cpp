#include <ROOT/RDataFrame.hxx>
#include <iostream>

#include "OversampledHistogram.h"

int main() {
  ROOT::RDataFrame frame(6);
  auto rdf = frame.Define("genEvent", "static_cast<unsigned long>(rdfentry_ / 2)").Define("value", "0.5");

  SequentialOversampledHistogram1D helper(2, "sequential", "Example", 2, 0, 2);
  auto histogram = rdf.Book<unsigned long, double>(helper, {"genEvent", "value"});
  std::cout << "Content: " << histogram->GetBinContent(1) << ", error: " << histogram->GetBinError(1) << '\n';
}
