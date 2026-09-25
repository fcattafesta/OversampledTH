#!/usr/bin/env python3
"""Measure exact and slot-local TH1 bin errors in one RDF event loop."""

import argparse
import json
from pathlib import Path

import ROOT


REPO_ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--factor", type=int, default=9)
    parser.add_argument("--spec", type=Path, default=REPO_ROOT / "data" / "site_sample.json")
    args = parser.parse_args()
    if args.threads < 1 or args.factor < 1:
        parser.error("threads and factor must be positive")

    spec = args.spec.resolve()
    sample = next(iter(json.loads(spec.read_text())["samples"].values()))
    for path in sample["files"]:
        if not Path(path).is_file():
            parser.error(f"Sample file is unavailable: {path}")

    if args.threads > 1:
        ROOT.EnableImplicitMT(args.threads)
    ROOT.gInterpreter.AddIncludePath(str(REPO_ROOT / "include"))
    ROOT.gInterpreter.Declare('#include "OversampledHistogram.h"')

    rdf = ROOT.RDF.Experimental.FromSpec(str(spec))
    rdf = (
        rdf.Define("nGenJet", "static_cast<unsigned int>(GenJet_pt.size())")
        .DefinePerSample("range_begin", "static_cast<ULong64_t>(rdfsampleinfo_.EntryRange().first)")
        .DefinePerSample("range_end", "static_cast<ULong64_t>(rdfsampleinfo_.EntryRange().second)")
        .Define("scaled_weight", f"1.0 / {args.factor}")
    )
    selected = rdf.Filter("ROOT::VecOps::Sum(Jet_pt > 30) >= 2")

    bookings = []
    for label, source, column in (
        ("nMuon", rdf, "nMuon"),
        ("nGenJet", rdf, "nGenJet"),
        ("selected_nMuon", selected, "nMuon"),
    ):
        exact_helper = ROOT.RangeAwareOversampledHistogram1D(
            args.factor, f"exact_{label}", label, 5, -0.5, 4.5
        )
        approx_helper = ROOT.SlotLocalOversampledHistogram1D(
            args.factor, f"approx_{label}", label, 5, -0.5, 4.5
        )
        exact = source.Book(exact_helper, ("event", column, "range_begin", "range_end"))
        approx = source.Book(approx_helper, ("event", column))
        rowwise = source.Histo1D((f"rowwise_{label}", label, 5, -0.5, 4.5), column, "scaled_weight")
        bookings.append((label, exact, approx, rowwise))

    ROOT.RDF.RunGraphs([result for _, exact, approx, rowwise in bookings for result in (exact, approx, rowwise)])

    rows = []
    for label, exact_result, approx_result, rowwise_result in bookings:
        exact = exact_result.GetValue()
        approx = approx_result.GetValue()
        rowwise = rowwise_result.GetValue()
        for bin_index in range(exact.GetNbinsX() + 2):
            content = exact.GetBinContent(bin_index)
            approx_content = approx.GetBinContent(bin_index)
            if abs(content - approx_content) > 1e-7 + 1e-10 * abs(content):
                raise RuntimeError(f"Content mismatch in {label} bin {bin_index}: {content} vs {approx_content}")
            rows.append(
                {
                    "histogram": label,
                    "bin": bin_index,
                    "bin_label": (
                        "underflow"
                        if bin_index == 0
                        else "overflow"
                        if bin_index == exact.GetNbinsX() + 1
                        else str(exact.GetBinCenter(bin_index))
                    ),
                    "content": content,
                    "exact_error": exact.GetBinError(bin_index),
                    "approx_error": approx.GetBinError(bin_index),
                    "rowwise_error": rowwise.GetBinError(bin_index),
                }
            )

    print("RESULT_JSON=" + json.dumps({"threads": args.threads, "bins": rows}, separators=(",", ":")))


if __name__ == "__main__":
    main()
