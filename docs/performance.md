# Performance and study results

## Reproducing the synthetic benchmark

Build `benchmark_actions`, then run:

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH="$(root-config --prefix)"
cmake --build build --target benchmark_actions
python studies/run_benchmarks.py --events 50000 --threads 1 2 4 --repeats 2
```

The runner writes raw measurements to [`studies/output/performance.csv`](../studies/output/performance.csv) and a plot to [`studies/output/performance.png`](../studies/output/performance.png) when Matplotlib is installed. Each fresh process processes 50,000 generated events, three rows per event, into one bin. The timed section includes row processing, worker startup, and finalization. It excludes action construction, ROOT file I/O, and RDF scheduling. `peak_rss_growth_kib` is a Linux `ru_maxrss` high-water increase after a warmup action; it is coarse and can be zero for small allocations.

## Observed rates

The checked-in run used ROOT 6.36.02 on 2026-09-25. Mean rates from two runs, rounded to millions of rows per second:

| Action | 1 thread | 2 threads | 4 threads |
| --- | ---: | ---: | ---: |
| Sequential | 6.49 | — | — |
| Slot-local | 1.74 | 1.49 | 1.30 |
| Range-aware | 2.87 | 1.35 | 1.19 |

This short synthetic workload did not benefit from more threads. The diagnostic ID sets, boundary merging, worker startup, and shared machine load all affect these rates. They are a comparison of these action paths on this host, not a prediction of full analysis throughput. Rerun with representative bins, weights, event sizes, and input files before choosing a deployment configuration.

## Exact and approximate errors by bin and thread count

`studies/compare_bin_errors.py` books three calculations on the **same RDF graph and event loop**: range-aware event grouping (exact under the ordering assumptions), slot-local event grouping (approximate when an event is fragmented), and conventional row-by-row weighted filling. It compares `nMuon`, `nGenJet`, and selected `nMuon`, including overflow. Each worker process uses one thread count; repeated processes expose scheduling variability.

```sh
python studies/compare_bin_errors.py --threads 1 2 4 8 16 32 --repeats 2
```

The checked-in run used the site sample, oversampling factor 9, ROOT 6.36.02, and two runs per thread count. The largest absolute **mean** slot-local versus range-aware bin-error difference among the three histograms was:

| Threads | Largest per-bin difference |
| ---: | ---: |
| 1 | 0% |
| 2 | 0.00126% |
| 4 | 0.00987% |
| 8 | 0.03034% |
| 16 | 0.06069% |
| 32 | 0.03034% |

The signed [per-bin plot](../studies/output/bin_errors.png), [raw measurements](../studies/output/bin_errors.csv), [per-bin summary](../studies/output/bin_errors_summary.csv), and [run diagnostics](../studies/output/bin_error_runs.csv) are in `studies/output/`. **Every observed individual bin difference was below 0.1%; the maximum was 0.06069%** in this two-run, 1–32-thread study. Most affected bins have a slightly smaller slot-local error. The trend is not strictly monotonic because RDF scheduling changes between runs. Some runs revisited an event within the same slot, so a zero cross-slot count alone is insufficient to establish exactness. This bound applies only to these three histograms and this sample; a sparse bin dominated by one fragmented event can differ more.

The separate [row-wise comparison](../studies/output/rowwise_vs_exact.png) shows the much larger uncertainty effect of **ignoring oversampling correlation itself**. For every `nGenJet` bin in this sample, the conventional row-wise error was about 66.7% smaller than the event-grouped error. Both custom actions apply event-level grouping; their much smaller gap measures the loss of correlation from fragmented processing, not the main oversampling correction.

## Slot-split study

`studies/scan_slot_splits.py` repeats the slot-local PyROOT example across thread counts. It writes [`slot_splits_scan.json`](../studies/output/slot_splits_scan.json) and [`slot_splits_scan.png`](../studies/output/slot_splits_scan.png) under `studies/output/`. The checked-in site-sample scan used five repetitions at 2, 5, 14, 36, 95, and 250 threads. Mean split counts were 1, 9.4, 31.2, 54.6, 67, and 67, respectively. These counts describe that sample and RDF scheduling; they are not universal limits.

```sh
python studies/scan_slot_splits.py --runs 5 --points 6 --min-threads 2 --max-threads 250
```
