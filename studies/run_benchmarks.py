#!/usr/bin/env python3
"""Measure action throughput on synthetic, event-aligned C++ workloads."""

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
import statistics
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build" / "benchmark_actions")
    parser.add_argument("--events", type=int, default=100_000)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, default=ROOT / "studies" / "output" / "performance.csv")
    args = parser.parse_args()
    if args.events < 1 or args.repeats < 1 or any(t < 1 for t in args.threads):
        parser.error("events, repeats, and thread counts must be positive")

    root_version = subprocess.check_output(["root-config", "--version"], text=True).strip()
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    rows = []
    for mode in ("sequential", "slot-local", "range-aware"):
        for threads in ([1] if mode == "sequential" else args.threads):
            for repeat in range(1, args.repeats + 1):
                proc = subprocess.run(
                    [str(args.binary), mode, str(args.events), str(threads)],
                    capture_output=True,
                    text=True,
                    check=True,
                )
                result = next((line for line in proc.stdout.splitlines() if line.startswith("RESULT,")), None)
                if result is None:
                    raise RuntimeError(f"Benchmark produced no RESULT row:\n{proc.stdout}\n{proc.stderr}")
                _, measured_mode, events, actual_threads, seconds, rows_per_second, rss_growth = result.split(",")
                row = {
                    "timestamp_utc": timestamp,
                    "root_version": root_version,
                    "mode": measured_mode,
                    "events": events,
                    "threads": actual_threads,
                    "repeat": repeat,
                    "seconds": seconds,
                    "rows_per_second": rows_per_second,
                    "peak_rss_growth_kib": rss_growth,
                }
                rows.append(row)
                print(f"{mode:11s} threads={threads} repeat={repeat}: {float(rows_per_second):,.0f} rows/s")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.output}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
    except ImportError:
        return

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for mode in ("sequential", "slot-local", "range-aware"):
        thread_counts = sorted({int(row["threads"]) for row in rows if row["mode"] == mode})
        rates = [
            statistics.mean(
                float(row["rows_per_second"])
                for row in rows
                if row["mode"] == mode and int(row["threads"]) == threads
            )
            for threads in thread_counts
        ]
        ax.plot(thread_counts, rates, marker="o", label=mode)
    ax.set(xlabel="Threads", ylabel="Rows per second", title=f"Synthetic action throughput ({args.events:,} events)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    plot = args.output.with_suffix(".png")
    fig.savefig(plot, dpi=160)
    print(f"Wrote {plot}")


if __name__ == "__main__":
    main()
