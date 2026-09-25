#!/usr/bin/env python3
"""Compare range-aware and slot-local bin errors across ROOT thread counts."""

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import statistics
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
WORKER = Path(__file__).with_name("measure_bin_errors.py")
SLOT_SPLITS_RE = re.compile(r"split across slots=(\d+)")
SLOT_REVISITS_RE = re.compile(r"revisited within slots=(\d+)")
RANGE_SPLITS_RE = re.compile(r"split across ranges=(\d+)")


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_worker(python, threads, factor, spec):
    proc = subprocess.run(
        [python, str(WORKER), "--threads", str(threads), "--factor", str(factor), "--spec", str(spec)],
        capture_output=True,
        text=True,
    )
    if proc.returncode:
        raise RuntimeError(f"Worker failed at {threads} threads:\n{proc.stdout}\n{proc.stderr}")
    payload = next((line.removeprefix("RESULT_JSON=") for line in proc.stdout.splitlines() if line.startswith("RESULT_JSON=")), None)
    if payload is None:
        raise RuntimeError(f"Worker produced no bin data:\n{proc.stdout}\n{proc.stderr}")
    return (
        json.loads(payload),
        [int(x) for x in SLOT_SPLITS_RE.findall(proc.stdout)],
        [int(x) for x in SLOT_REVISITS_RE.findall(proc.stdout)],
        [int(x) for x in RANGE_SPLITS_RE.findall(proc.stdout)],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8, 16, 32])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--factor", type=int, default=9)
    parser.add_argument("--spec", type=Path, default=REPO_ROOT / "data" / "site_sample.json")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "studies" / "output")
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()
    if args.repeats < 1 or args.factor < 1 or not args.threads or any(t < 1 for t in args.threads):
        parser.error("threads, repeats, and factor must be positive")
    threads = sorted(set(args.threads))
    spec = args.spec.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    root_version = subprocess.check_output(["root-config", "--version"], text=True).strip()
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    raw = []
    runs = []
    total = len(threads) * args.repeats
    for index, thread_count in enumerate(threads):
        for repeat in range(1, args.repeats + 1):
            print(f"[{index * args.repeats + repeat}/{total}] threads={thread_count}, repeat={repeat}", flush=True)
            result, slot_splits, slot_revisits, range_splits = run_worker(
                args.python, thread_count, args.factor, spec
            )
            runs.append(
                {
                    "timestamp_utc": timestamp,
                    "root_version": root_version,
                    "threads": thread_count,
                    "repeat": repeat,
                    "max_cross_slot_event_ids": max(slot_splits, default=0),
                    "max_same_slot_revisited_event_ids": max(slot_revisits, default=0),
                    "max_multi_range_event_ids": max(range_splits, default=0),
                }
            )
            for bin_row in result["bins"]:
                exact = bin_row["exact_error"]
                approx = bin_row["approx_error"]
                rowwise = bin_row["rowwise_error"]
                raw.append(
                    {
                        "threads": thread_count,
                        "repeat": repeat,
                        "histogram": bin_row["histogram"],
                        "bin": bin_row["bin"],
                        "bin_label": bin_row["bin_label"],
                        "content": bin_row["content"],
                        "exact_error": exact,
                        "approx_error": approx,
                        "rowwise_error": rowwise,
                        "error_difference": approx - exact,
                        "relative_difference_pct": 100 * (approx / exact - 1) if exact else "",
                        "rowwise_relative_difference_pct": 100 * (rowwise / exact - 1) if exact else "",
                    }
                )
            print(
                f"    cross-slot IDs={runs[-1]['max_cross_slot_event_ids']}, "
                f"same-slot revisits={runs[-1]['max_same_slot_revisited_event_ids']}, "
                f"multi-range IDs={runs[-1]['max_multi_range_event_ids']}",
                flush=True,
            )

    groups = {}
    for row in raw:
        key = row["threads"], row["histogram"], row["bin"], row["bin_label"]
        groups.setdefault(key, []).append(row)
    summary = []
    for (thread_count, histogram, bin_index, bin_label), values in sorted(groups.items()):
        relative = [float(row["relative_difference_pct"]) for row in values if row["relative_difference_pct"] != ""]
        rowwise_relative = [
            float(row["rowwise_relative_difference_pct"])
            for row in values
            if row["rowwise_relative_difference_pct"] != ""
        ]
        summary.append(
            {
                "threads": thread_count,
                "histogram": histogram,
                "bin": bin_index,
                "bin_label": bin_label,
                "runs": len(values),
                "mean_content": statistics.mean(float(row["content"]) for row in values),
                "mean_exact_error": statistics.mean(float(row["exact_error"]) for row in values),
                "mean_approx_error": statistics.mean(float(row["approx_error"]) for row in values),
                "mean_rowwise_error": statistics.mean(float(row["rowwise_error"]) for row in values),
                "mean_relative_difference_pct": statistics.mean(relative) if relative else "",
                "min_relative_difference_pct": min(relative) if relative else "",
                "max_relative_difference_pct": max(relative) if relative else "",
                "mean_rowwise_relative_difference_pct": (
                    statistics.mean(rowwise_relative) if rowwise_relative else ""
                ),
            }
        )

    write_csv(out_dir / "bin_errors.csv", raw)
    write_csv(out_dir / "bin_errors_summary.csv", summary)
    write_csv(out_dir / "bin_error_runs.csv", runs)

    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    histograms = ("nMuon", "nGenJet", "selected_nMuon")
    fig, axes = plt.subplots(len(histograms), 1, figsize=(9, 10), sharex=True, constrained_layout=True)
    for ax, histogram in zip(axes, histograms):
        bins = sorted({int(row["bin"]) for row in summary if row["histogram"] == histogram})
        for bin_index in bins:
            points = [
                row for row in summary
                if row["histogram"] == histogram and int(row["bin"]) == bin_index
                and row["mean_relative_difference_pct"] != ""
            ]
            if not points:
                continue
            if all(abs(float(row["mean_relative_difference_pct"])) < 1e-10 for row in points):
                continue
            ax.plot(
                [int(row["threads"]) for row in points],
                [float(row["mean_relative_difference_pct"]) for row in points],
                marker="o",
                label=f"bin {points[0]['bin_label']}",
            )
        ax.axhline(0, color="black", linewidth=0.7)
        ax.set_title(histogram)
        ax.set_ylabel("Approx − exact error (%)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, ncol=3)
    axes[-1].set_xticks(threads, [str(value) for value in threads])
    axes[-1].set_xlabel("ROOT threads")
    fig.suptitle(f"Per-bin error difference: slot-local vs range-aware ({args.repeats} runs/point)")
    fig.savefig(out_dir / "bin_errors.png", dpi=160)

    fig, axes = plt.subplots(len(histograms), 1, figsize=(9, 9), constrained_layout=True)
    baseline_thread = threads[0]
    for ax, histogram in zip(axes, histograms):
        points = [
            row for row in summary
            if row["histogram"] == histogram and int(row["threads"]) == baseline_thread
            and row["mean_rowwise_relative_difference_pct"] != ""
        ]
        ax.bar(
            [row["bin_label"] for row in points],
            [float(row["mean_rowwise_relative_difference_pct"]) for row in points],
            color="#4C78A8",
        )
        ax.axhline(0, color="black", linewidth=0.7)
        ax.set_title(histogram)
        ax.set_ylabel("Row-wise − exact error (%)")
        ax.grid(axis="y", alpha=0.25)
    axes[-1].set_xlabel("Histogram bin")
    fig.suptitle("Oversampling correlation effect: independent rows vs grouped events")
    fig.savefig(out_dir / "rowwise_vs_exact.png", dpi=160)

    print(f"Wrote bin-error CSVs and plot to {out_dir}")
    for thread_count in threads:
        relevant = [
            abs(float(row["mean_relative_difference_pct"]))
            for row in summary
            if int(row["threads"]) == thread_count and row["mean_relative_difference_pct"] != ""
        ]
        print(f"threads={thread_count}: largest |mean per-bin error difference|={max(relevant, default=0):.6g}%")


if __name__ == "__main__":
    main()
