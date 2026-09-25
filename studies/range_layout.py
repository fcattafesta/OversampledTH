#!/usr/bin/env python3
"""Inspect the RDF entry ranges that determine boundary-event spills."""

import argparse
from pathlib import Path

import ROOT

DEFAULT_SPEC = Path(__file__).resolve().parents[1] / "data" / "site_sample.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-j", type=int, default=1, help="Number of ROOT threads")
    parser.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    args = parser.parse_args()

    if args.j > 1:
        ROOT.EnableImplicitMT(args.j)

    rdf = ROOT.RDF.Experimental.FromSpec(str(args.spec.resolve()))
    rdf = rdf.DefinePerSample(
        "range_key",
        "std::to_string(static_cast<ULong64_t>(rdfsampleinfo_.EntryRange().first)) + ':' + "
        "std::to_string(static_cast<ULong64_t>(rdfsampleinfo_.EntryRange().second))",
    )
    keys = rdf.Take["std::string"]("range_key").GetValue()
    ranges = sorted({tuple(map(int, str(key).split(":"))) for key in keys})

    print(" idx | begin    | end      | span")
    print("-----+----------+----------+---------")
    for i, (begin, end) in enumerate(ranges):
        print(f" {i:>3} | {begin:>8} | {end:>8} | {end - begin:>7}")

    overlaps = sum(b1 < e0 for (_, e0), (b1, _) in zip(ranges, ranges[1:]))
    gaps = sum(b1 > e0 for (_, e0), (b1, _) in zip(ranges, ranges[1:]))
    print(f"Ranges: {len(ranges)}; overlaps/backwards: {overlaps}; gaps: {gaps}")


if __name__ == "__main__":
    main()
