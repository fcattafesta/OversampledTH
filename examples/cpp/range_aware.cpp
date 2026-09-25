#include <iostream>

#include "OversampledHistogram.h"

int main() {
  ROOT::EnableImplicitMT(2);
  RangeAwareOversampledHistogram1D helper(2, "range_aware", "Example", 2, 0, 2);

  // One generated event is split between two entry ranges and two slots.
  // In an RDF job, pass these range columns from rdfsampleinfo_.EntryRange().
  helper.Exec(0, 7UL, 0.5, 0UL, 2UL, 2.0);
  helper.Exec(1, 7UL, 0.5, 2UL, 4UL, 4.0);
  helper.Finalize();

  const auto histogram = helper.GetResultPtr();
  std::cout << "Content: " << histogram->GetBinContent(1) << ", error: " << histogram->GetBinError(1) << '\n';
}
