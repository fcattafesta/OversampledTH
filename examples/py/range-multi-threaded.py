import argparse
import ROOT
from common import DEFAULT_SPEC, REPO_ROOT, sample_from_spec

ROOT.gInterpreter.AddIncludePath(str(REPO_ROOT / "include"))
ROOT.gInterpreter.Declare('#include "OversampledHistogram.h"')

parser = argparse.ArgumentParser(description="Range-aware oversampled histogram on a NanoAOD sample")
parser.add_argument("-j", type=int, default=1, help="Number of threads to use")
parser.add_argument("--spec", default=DEFAULT_SPEC, help="RDF sample spec JSON")
args = parser.parse_args()
sample_from_spec(args.spec)

if args.j > 1:
    print(f"Starting range test with implicit multi-threading enabled. Using {args.j} threads.")
    ROOT.EnableImplicitMT(args.j)
else:
    print("Starting range test with implicit multi-threading disabled.")

rdf = ROOT.RDF.Experimental.FromSpec(str(args.spec))
rdf = (
    rdf.Define("nGenJet", "static_cast<unsigned int>(GenJet_pt.size())")
    .DefinePerSample("range_begin", "static_cast<ULong64_t>(rdfsampleinfo_.EntryRange().first)")
    .DefinePerSample("range_end", "static_cast<ULong64_t>(rdfsampleinfo_.EntryRange().second)")
)
rdf_base = rdf.Filter("fold == 0")
rdf_sel = rdf.Filter("ROOT::VecOps::Sum(Jet_pt > 30) >= 2")
rdf_base_sel = rdf_base.Filter("ROOT::VecOps::Sum(Jet_pt > 30) >= 2")

h_base = rdf_base.Histo1D(("h_base_nMuon", "Base nMuon Histogram;Number of Muons;Entries", 5, -0.5, 4.5), "nMuon")

helper_0 = ROOT.RangeAwareOversampledHistogram1D(9, "h_oversampled", "Range Oversampled nMuon Histogram", 5, -0.5, 4.5)
h_oversampled = rdf.Book(helper_0, ("event", "nMuon", "range_begin", "range_end"))

helper_1 = ROOT.RangeAwareOversampledHistogram1D(9, "h_gen_over", "Range Oversampled Gen Histogram", 5, -0.5, 4.5)
h_gen_over = rdf.Book(helper_1, ("event", "nGenJet", "range_begin", "range_end"))

helper_2 = ROOT.RangeAwareOversampledHistogram1D(9, "h_fold1_over", "Range Oversampled nMuon (nJet >= 30) >= 2", 5, -0.5, 4.5)
h_fold1_over = rdf_sel.Book(helper_2, ("event", "nMuon", "range_begin", "range_end"))

h_gen_base = rdf_base.Histo1D(("h_gen_base_1", "Base Gen nJet Histogram;Number of Gen Jets;Entries", 5, -0.5, 4.5), "nGenJet")
h_fold1_base = rdf_base_sel.Histo1D(("h_fold1_base", "Base nMuon (nJet >= 30) >= 2;Number of Muons;Entries", 5, -0.5, 4.5), "nMuon")

print("== nMuon ==")
h_oversampled.Print("all")
h_base.Sumw2()
h_base.Print("all")

print("\n== nGenJet ==")
h_gen_over.Print("all")
h_gen_base.Sumw2()
h_gen_base.Print("all")

print("\n== nMuon (nJet >= 30) >= 2 ==")
h_fold1_over.Print("all")
h_fold1_base.Sumw2()
h_fold1_base.Print("all")
