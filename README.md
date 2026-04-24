# OversampledTH

A ROOT RDataFrame action for creating histograms with proper uncertainty handling when performing oversampling in Flashsim.

This repository currently includes:

- `STOversampledTH`: single-thread helper
- `DumbOversampledTH`: multi-thread helper that assumes the same `genEvent` is not split across slots

The dumb MT helper prints diagnostics in `Finalize()`:

- number of unique `genEvent`
- number of `genEvent` seen in more than one slot (`split across slots`)

The split count quantifies the amount of assumption violation in MT execution.

## Testing

Run the test suite to compare performance and accuracy:

### Single-thread test

```sh
python test/single-threaded.py
```

### Dumb multi-thread test

```sh
python test/dumb-multi-threaded.py -j 4
```

This prints histogram summaries and the dumb MT diagnostics line:

```text
DumbOversampledTH diagnostics: unique genEvents=..., split across slots=...
```

### Scan split-fault distribution vs threads

Use the scan utility to measure split-fault behavior as a function of thread count, with repeated runs for uncertainty estimation.

Default configuration:

- 6 log-spaced thread points from 2 to 250
- 5 repetitions per point

Run:

```sh
python test/scan_dumb_mt_faults.py
```

Optional overrides:

```sh
python test/scan_dumb_mt_faults.py --runs 5 --points 6 --min-threads 2 --max-threads 250
```

Outputs:

- `test/dumb_mt_faults_scan.png`: plot with mean/std vs thread count and per-thread distributions
- `test/dumb_mt_faults_scan.json`: raw samples, means, and standard deviations