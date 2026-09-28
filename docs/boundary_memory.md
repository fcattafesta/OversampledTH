# Sparse range-boundary storage

`RangeAwareOversampledHistogram` stores deferred events as sorted vectors of nonzero bin contents instead of ROOT histogram clones. Each record includes its boundary-piece count. Finalization moves the first record for an event into the global map, merges subsequent records by bin index, and releases consumed slot entries and map buckets. No dense global boundary histogram is created.

## Uncertainty semantics

All pieces of one event are combined before the output histogram receives one native `TH1::Fill` per nonzero bin, with weight `total / factor`. ROOT `Sumw2` therefore retains the event-level variance. Intermediate histogram error arrays were never used by this calculation and are not needed in sparse storage. Underflow and overflow bins are included; exact zero totals from cancellation are omitted.

Sparse contents use `float` for `TH1F` and `double` for `TH1D`, retaining the supported histogram types' intermediate precision. Numerical equivalence is tested rather than universal bitwise equivalence: event summation order and ROOT floating-point behavior can affect the last digits.

## Memory and time limits

For `S` slots, `B` bins, and `K` occupied bins retained across deferred event records, storage is approximately `O(S * B + K)` plus map metadata. Each slot still owns two dense working histograms. Accumulator merging temporarily allocates a vector proportional to the two sparse inputs, rather than another full histogram.

The current event is still accumulated in a ROOT histogram. Extracting a boundary contribution scans its bins, so this change does not remove the `O(B)` extraction cost. Deferred storage can still grow with the number of ranges and events; it is not a bounded working set. At full occupancy, sparse records have little memory advantage and can have additional container overhead.

## Measured comparison

The [CSV](../studies/output/boundary_memory.csv) and [plot](../studies/output/boundary_memory.png) compare the dense implementation at commit `59dc246` with sparse storage using ROOT 6.36.02 and GCC 13, optimized builds, on 2026-09-28. Each case runs in a fresh process after ROOT warmup. RSS columns are Linux `ru_maxrss` high-water measurements at construction, processing, and finalization; the reported growth includes working histograms and merge-time allocations. Zero growth means the workload did not exceed the previous high-water mark, not that it allocated no memory.

Selected peak-RSS growth results:

| Events | Factor | Bins | Slots | Range rows | Occupied bins/event | Dense MiB | Sparse MiB |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 9 | 10,000 | 8 | 7 | 1 | 372 | 2 |
| 1,000 | 9 | 10,000 | 8 | 7 | 20 | 374 | 4 |
| 1,000 | 9 | 10,000 | 8 | 257 | 1 | 16 | 4 |
| 100 | 1,000 | 10,000 | 32 | 257 | 1 | 88 | 10 |
| 200 | 9 | 1,000 | 8 | 7 | 1,000 | 8 | 8 |

These are controlled serialized action calls with deterministic round-robin range assignment, not an RDF scheduling or parallel-throughput benchmark. All occupied output bins must have the expected content, event-level error, and fill count before a measurement is accepted.

## Reproduction

With a ROOT environment loaded:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$(root-config --prefix)"
cmake --build build --target benchmark_boundary_memory
python studies/run_boundary_memory.py --binary build/benchmark_boundary_memory
```

To compare with the old dense header without changing the working checkout:

```sh
mkdir -p /tmp/oversampled-dense/include
git show 59dc246:include/OversampledHistogram.h > /tmp/oversampled-dense/include/OversampledHistogram.h
c++ -O3 -std=c++17 -I/tmp/oversampled-dense/include studies/benchmark_boundary_memory.cpp \
    $(root-config --cflags --libs) -o /tmp/oversampled-dense/benchmark_boundary_memory
python studies/run_boundary_memory.py --binary build/benchmark_boundary_memory \
    --baseline-binary /tmp/oversampled-dense/benchmark_boundary_memory
```

`test_sparse_boundaries.cpp` checks dense-reference equivalence for both supported precisions and enforces a 64 MiB peak-RSS growth budget for eight slots, 10,000 bins, and 1,000 generated events. The existing concurrent coherence tests continue to exercise split events across worker threads.
