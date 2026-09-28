#!/usr/bin/env python3
"""Compare boundary-memory growth with controlled bin occupancy and range assignment."""

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
import subprocess

REPO_ROOT = Path(__file__).resolve().parents[1]
FIELDS = [
    "events", "factor", "bins", "slots", "range_rows", "occupancy",
    "fill_seconds", "finalize_seconds", "rss_before_kib", "rss_slots_kib", "rss_filled_kib", "peak_rss_kib",
]
# Independent changes to bins, slots, occupancy, factor, and range size.
CASES = [
    (1000, 9, 100, 8, 7, 1),
    (1000, 9, 10000, 1, 7, 1),
    (1000, 9, 10000, 8, 7, 1),
    (1000, 9, 10000, 8, 7, 20),
    (1000, 9, 10000, 8, 257, 1),
    (100, 1000, 10000, 8, 257, 1),
    (100, 1000, 10000, 32, 257, 1),
    (200, 9, 1000, 8, 7, 1000),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=REPO_ROOT / "build" / "benchmark_boundary_memory")
    parser.add_argument("--baseline-binary", type=Path, help="Optional executable built against the old dense header")
    parser.add_argument("--output", type=Path, default=REPO_ROOT / "studies" / "output" / "boundary_memory.csv")
    args = parser.parse_args()
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    version = subprocess.check_output(["root-config", "--version"], text=True).strip()
    binaries = [("sparse", args.binary)]
    if args.baseline_binary:
        binaries.insert(0, ("dense", args.baseline_binary))
    results = []
    for implementation, binary in binaries:
        for case in CASES:
            proc = subprocess.run([str(binary.resolve()), *map(str, case)], capture_output=True, text=True, timeout=120)
            if proc.returncode:
                raise RuntimeError(f"{implementation} {case} failed:\n{proc.stdout}\n{proc.stderr}")
            line = next((line for line in proc.stdout.splitlines() if line.startswith("RESULT,")), None)
            if line is None:
                raise RuntimeError(f"Missing RESULT row: {proc.stdout}")
            row = {"timestamp_utc": timestamp, "root_version": version, "implementation": implementation}
            row.update(zip(FIELDS, line.split(",")[1:]))
            row["peak_rss_growth_kib"] = int(row["peak_rss_kib"]) - int(row["rss_before_kib"])
            results.append(row)
            print(f"{implementation} {case}: {row['peak_rss_growth_kib'] / 1024:.2f} MiB peak RSS growth", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    print(f"Wrote {args.output}")
    plot_results(args.output, results, [name for name, _ in binaries])


def plot_results(output, results, implementations):
    try:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)
    width = 0.8 / len(implementations)
    for offset, implementation in enumerate(implementations):
        values = [float(row["peak_rss_growth_kib"]) / 1024 for row in results if row["implementation"] == implementation]
        positions = [i - 0.4 + width * (offset + 0.5) for i in range(len(CASES))]
        bars = ax.bar(positions, values, width=width, label=implementation)
        ax.bar_label(bars, fmt="%.1f", fontsize=8)
    labels = [f"E={e}, F={f}\nB={b}, S={s}\nR={r}, K={o}" for e, f, b, s, r, o in CASES]
    ax.set_xticks(range(len(CASES)), labels, fontsize=8)
    ax.set(ylabel="Peak RSS growth after warmup (MiB)", title="Deferred boundary storage: dense vs sparse",
           xlabel="E: events; F: oversampling factor; B: bins; S: slots; R: rows per range; K: occupied bins per event")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.savefig(output.with_suffix(".png"), dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
