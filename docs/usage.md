# Using the histogram actions

Include `OversampledHistogram.h` and book one of the `1D` or `1F` aliases. The constructor is `(factor, name, title, number_of_bins, xmin, xmax)`; the factor must be at least one. The first booked column is the generated-event ID, the second is a scalar or `ROOT::VecOps::RVec` value, and an optional last column is the row weight.

## Sequential RDF

Keep implicit MT disabled and ensure equal generated-event IDs are contiguous:

```python
import ROOT

ROOT.gInterpreter.AddIncludePath("include")
ROOT.gInterpreter.Declare('#include "OversampledHistogram.h"')
helper = ROOT.SequentialOversampledHistogram1D(9, "h", "nMuon", 5, -0.5, 4.5)
histogram = rdf.Book(helper, ("event", "nMuon"))
histogram.GetValue()  # Runs the RDF event loop.
```

The complete Python example is [`examples/py/single-threaded.py`](../examples/py/single-threaded.py); the C++ version is [`examples/cpp/sequential.cpp`](../examples/cpp/sequential.cpp).

## Slot-local RDF

Call `ROOT.EnableImplicitMT(n)` before constructing the helper and use `SlotLocalOversampledHistogram1D`. Its `Finalize()` diagnostic counts generated-event IDs seen in more than one slot **and** IDs revisited after another event in the same slot. Either count can indicate that the histogram content is scaled but its bin errors omit correlations between fragments. See the [Python](../examples/py/slot-local-multi-threaded.py) and [C++](../examples/cpp/slot_local.cpp) examples.

## Range-aware RDF

Use a sample-spec RDF so `rdfsampleinfo_.EntryRange()` is available. Define the bounds per sample and pass them after the value column:

```python
rdf = ROOT.RDF.Experimental.FromSpec("data/site_sample.json")
rdf = (
    rdf.DefinePerSample("range_begin", "static_cast<ULong64_t>(rdfsampleinfo_.EntryRange().first)")
       .DefinePerSample("range_end", "static_cast<ULong64_t>(rdfsampleinfo_.EntryRange().second)")
)
helper = ROOT.RangeAwareOversampledHistogram1D(9, "h", "nMuon", 5, -0.5, 4.5)
histogram = rdf.Book(helper, ("event", "nMuon", "range_begin", "range_end"))
```

Enable implicit MT before creating `rdf` and `helper` when using multiple threads. The full [Python example](../examples/py/range-multi-threaded.py) books an RDF; the [C++ example](../examples/cpp/range_aware.cpp) demonstrates a controlled event split across two ranges and slots.

## Inputs and commands

The Python examples default to [`data/site_sample.json`](../data/site_sample.json), which points to a site-local ROOT file. They accept `--spec` to use another RDF sample spec. Run from any working directory:

```sh
python examples/py/single-threaded.py --spec data/site_sample.json
python examples/py/slot-local-multi-threaded.py -j 4 --spec data/site_sample.json
python examples/py/range-multi-threaded.py -j 4 --spec data/site_sample.json
```

After building, the corresponding C++ executables are `build/example_sequential`, `build/example_slot_local`, and `build/example_range_aware`.
