# Studies

- `run_benchmarks.py` runs `benchmark_actions` in fresh processes and writes raw CSV and a throughput plot to `output/`.
- `run_boundary_memory.py` runs `benchmark_boundary_memory` across controlled bin, slot, factor, range-size, and occupancy cases; it optionally compares an old dense executable and writes CSV/PNG results.
- `compare_bin_errors.py` books range-aware, slot-local, and row-wise histograms in the same RDF graph at each thread count. It writes per-bin CSVs and plots to `output/`.
- `scan_slot_splits.py` measures slot-local split counts on the sample in `data/site_sample.json`; it writes JSON and a plot to `output/`.
- `range_layout.py` lists RDF entry ranges for the sample.

```sh
python studies/run_benchmarks.py --binary build/benchmark_actions
python studies/run_boundary_memory.py --binary build/benchmark_boundary_memory
python studies/compare_bin_errors.py --threads 1 2 4 8 16 32 --repeats 2
python studies/scan_slot_splits.py --runs 5 --points 6
python studies/range_layout.py -j 4
```

The throughput and boundary-memory benchmarks have no file input. The other commands require the site-local sample or a replacement `--spec`. See [`docs/performance.md`](../docs/performance.md) and [`docs/boundary_memory.md`](../docs/boundary_memory.md) for measured results and limits.
