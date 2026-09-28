# Automated tests

`test_histograms.cpp` checks weighted and vector fills, errors, underflow/overflow, copies, invalid factors, RDF booking, controlled slot/range splits, and same-slot revisit detection. `test_coherence.cpp` compares all actions with an independent event-level reference histogram and runs an event split across two real worker threads. `test_stress.cpp` processes 50,000 generated events and checks contents, errors, entry counts, and a Linux peak-RSS growth budget.

Run all tests through CTest after building:

```sh
ctest --test-dir build --output-on-failure
```

These tests use synthetic data and do not depend on the site-local file in `data/`.

`test_sparse_boundaries.cpp` compares sparse boundary merging against the former dense ROOT-histogram calculation for `TH1D` and `TH1F`, including cancellation, underflow/overflow, empty vectors, full occupancy, and repeated pieces across slots and ranges. It also processes 1,000 events with eight slots and 10,000 bins, enforcing a 64 MiB peak-RSS growth budget through finalization.
