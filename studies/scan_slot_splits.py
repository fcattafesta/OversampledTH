#!/usr/bin/env python3
import argparse
import json
import math
import os
import re
import statistics
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt


SPLIT_RE = re.compile(r"split across slots=(\d+)")


def build_logspace_threads(start: int, stop: int, n_points: int) -> list[int]:
    if n_points < 2:
        return [max(1, int(start))]

    vals: list[int] = []
    for i in range(n_points):
        t = i / (n_points - 1)
        v = math.exp(math.log(start) + t * (math.log(stop) - math.log(start)))
        vals.append(int(round(v)))

    vals = sorted(set(vals))
    if vals[0] != start:
        vals[0] = start
    if vals[-1] != stop:
        vals[-1] = stop

    while len(vals) < n_points:
        for i in range(1, len(vals)):
            mid = (vals[i - 1] + vals[i]) // 2
            if mid not in vals and mid >= start and mid <= stop:
                vals.append(mid)
                break
        vals = sorted(set(vals))

    if len(vals) > n_points:
        idxs = [round(i * (len(vals) - 1) / (n_points - 1)) for i in range(n_points)]
        vals = [vals[i] for i in idxs]

    return vals


def parse_split_count(output: str) -> int:
    matches = [int(x) for x in SPLIT_RE.findall(output)]
    if not matches:
        raise RuntimeError("Could not find 'split across slots=...' in command output")

    if len(set(matches)) > 1:
        print(
            f"Warning: found multiple split-count values in one run: {matches}; using max={max(matches)}",
            file=sys.stderr,
        )
    return max(matches)


def run_one(test_script: str, n_threads: int, python_bin: str, spec: str | None) -> tuple[int, str]:
    cmd = [python_bin, test_script, "-j", str(n_threads)]
    if spec:
        cmd.extend(["--spec", spec])
    proc = subprocess.run(cmd, capture_output=True, text=True)
    merged = (proc.stdout or "") + "\n" + (proc.stderr or "")

    if proc.returncode != 0:
        raise RuntimeError(
            f"Run failed for -j {n_threads} with exit code {proc.returncode}\n"
            f"Command: {' '.join(cmd)}\n"
            f"Output:\n{merged}"
        )

    return parse_split_count(merged), merged


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Scan SlotLocalOversampledHistogram split counts vs number of threads and plot summary."
    )
    parser.add_argument("--runs", type=int, default=5, help="Number of repetitions per thread count (default: 5)")
    parser.add_argument("--points", type=int, default=6, help="Number of logspace thread points (default: 6)")
    parser.add_argument("--min-threads", type=int, default=2, help="Minimum number of threads (default: 2)")
    parser.add_argument("--max-threads", type=int, default=250, help="Maximum number of threads (default: 250)")
    parser.add_argument(
        "--python-bin", default=sys.executable, help="Python executable used to run slot-local-multi-threaded.py"
    )
    parser.add_argument(
        "--test-script",
        default=str(Path(__file__).resolve().parents[1] / "examples" / "py" / "slot-local-multi-threaded.py"),
        help="Path to the slot-local example",
    )
    parser.add_argument("--spec", help="Optional RDF sample spec passed to the example")
    parser.add_argument(
        "--out-png",
        default=str(Path(__file__).parent / "output" / "slot_splits_scan.png"),
        help="Output PNG path",
    )
    parser.add_argument(
        "--out-json",
        default=str(Path(__file__).parent / "output" / "slot_splits_scan.json"),
        help="Output JSON path",
    )
    args = parser.parse_args()

    if args.min_threads < 1 or args.max_threads < args.min_threads:
        raise ValueError("Invalid thread range")
    if args.runs < 1:
        raise ValueError("--runs must be >= 1")
    if args.points < 1:
        raise ValueError("--points must be >= 1")
    if args.points > args.max_threads - args.min_threads + 1:
        raise ValueError("--points exceeds the available integer thread counts")

    thread_values = build_logspace_threads(args.min_threads, args.max_threads, args.points)
    print(f"Thread points: {thread_values}")

    samples: dict[int, list[int]] = {t: [] for t in thread_values}

    total_runs = len(thread_values) * args.runs
    run_idx = 0
    for t in thread_values:
        for r in range(args.runs):
            run_idx += 1
            print(f"[{run_idx}/{total_runs}] Running -j {t}, repetition {r + 1}/{args.runs} ...")
            splits, _ = run_one(args.test_script, t, args.python_bin, args.spec)
            samples[t].append(splits)
            print(f"    splits={splits}")

    means = [statistics.mean(samples[t]) for t in thread_values]
    stds = [statistics.pstdev(samples[t]) for t in thread_values]

    results = {
        "thread_values": thread_values,
        "runs_per_thread": args.runs,
        "samples": {str(k): v for k, v in samples.items()},
        "mean_splits": {str(t): means[i] for i, t in enumerate(thread_values)},
        "std_splits": {str(t): stds[i] for i, t in enumerate(thread_values)},
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out_json)), exist_ok=True)
    with open(args.out_json, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13, 5.5))

    x = list(range(len(thread_values)))
    ax0.bar(x, means, yerr=stds, capsize=6, color="#4C78A8", alpha=0.85)
    for i, t in enumerate(thread_values):
        ys = samples[t]
        xs = [i + (j - (len(ys) - 1) / 2) * 0.04 for j in range(len(ys))]
        ax0.plot(xs, ys, "o", color="#1f1f1f", markersize=4, alpha=0.75)
    ax0.set_xticks(x)
    ax0.set_xticklabels([str(t) for t in thread_values])
    ax0.set_xlabel("Threads")
    ax0.set_ylabel("Split genEvents")
    ax0.set_title("Split count vs threads (mean +/- std)")
    ax0.grid(axis="y", alpha=0.3)

    max_splits = max(max(v) for v in samples.values())
    bins = list(range(0, max_splits + 2))
    hist_data = [samples[t] for t in thread_values]
    labels = [f"j={t}" for t in thread_values]
    ax1.hist(hist_data, bins=bins, histtype="step", linewidth=1.8, label=labels)
    ax1.set_xlabel("Split genEvents")
    ax1.set_ylabel("Count across repeated runs")
    ax1.set_title("Split-count distributions")
    ax1.grid(axis="y", alpha=0.3)
    ax1.legend(fontsize=8)

    fig.suptitle(
        f"SlotLocalOversampledHistogram split scan ({args.runs} runs each, {len(thread_values)} thread points)",
        fontsize=12,
    )
    fig.tight_layout()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_png)), exist_ok=True)
    fig.savefig(args.out_png, dpi=160)
    print(f"Saved plot: {args.out_png}")
    print(f"Saved data: {args.out_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
