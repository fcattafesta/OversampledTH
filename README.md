# OversampledTH

Header-only C++17 histogram actions for ROOT `RDataFrame::Book`. They combine oversampled rows by generated event before filling a `TH1`, so bin errors reflect event-level contributions.

| Action | Execution | Required grouping |
| --- | --- | --- |
| `SequentialOversampledHistogram` | One thread | Equal event IDs are contiguous. |
| `SlotLocalOversampledHistogram` | Multiple RDF slots | Equal IDs are contiguous within one slot and never span slots. |
| `RangeAwareOversampledHistogram` | Multiple RDF slots | Equal IDs are contiguous in global entry order; pass RDF entry-range bounds. |

The public header is [`include/OversampledHistogram.h`](include/OversampledHistogram.h). The repository is organized as follows:

| Directory | Contents |
| --- | --- |
| [`examples/`](examples/) | PyROOT and C++ examples for each action. |
| [`data/`](data/) | Input sample specification. |
| [`test/`](test/) | Automated correctness, coherence, event-count, and memory tests. |
| [`studies/`](studies/) | Throughput benchmark, per-bin exact/approximate error comparison, slot-split scan, range inspection, and generated results in `output/`. |
| [`docs/`](docs/) | [Usage](docs/usage.md), [assumptions and practical tips](docs/assumptions.md), [performance](docs/performance.md), and [architecture](docs/architecture.md). |

## Build and test

ROOT with `Hist` and `ROOTDataFrame`, CMake, and a C++17 compiler are required. The automated tests use synthetic data:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH="$(root-config --prefix)"
cmake --build build
ctest --test-dir build --output-on-failure
```

The sample-specific PyROOT examples use the file listed in `data/site_sample.json`. Pass `--spec path/to/spec.json` to use another sample. Benchmark instructions and the scope of the measured results are in [`docs/performance.md`](docs/performance.md).
