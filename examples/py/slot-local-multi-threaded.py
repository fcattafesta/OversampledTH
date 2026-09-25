import argparse
import ROOT
from common import DEFAULT_SPEC, REPO_ROOT, sample_from_spec

ROOT.gInterpreter.AddIncludePath(str(REPO_ROOT / "include"))
ROOT.gInterpreter.Declare('#include "OversampledHistogram.h"')

parser = argparse.ArgumentParser(description="Slot-local oversampled histogram on a NanoAOD sample")
parser.add_argument("-j", type=int, default=1, help="Number of threads to use (default: 1)")
parser.add_argument("--spec", default=DEFAULT_SPEC, help="RDF sample spec JSON")
args = parser.parse_args()

if args.j > 1:
    print(f"Starting test with implicit multi-threading enabled. Using {args.j} threads.")
    ROOT.EnableImplicitMT(args.j)
else:
    print("Starting test with implicit multi-threading disabled.")

root_file, tree_name = sample_from_spec(args.spec)
file = ROOT.TFile.Open(root_file)
events = file.Get(tree_name)
print(f"Number of events: {events.GetEntries()}")

rdf = ROOT.RDataFrame(tree_name, root_file).Define("nGenJet", "static_cast<unsigned int>(GenJet_pt.size())")
rdf_base = rdf.Filter("fold == 0")
rdf_sel = rdf.Filter("ROOT::VecOps::Sum(Jet_pt > 30) >= 2")
rdf_base_sel = rdf_base.Filter("ROOT::VecOps::Sum(Jet_pt > 30) >= 2")

h_base = rdf_base.Histo1D(("h_base_nMuon", "Base nMuon Histogram;Number of Muons;Entries", 5, -0.5, 4.5), "nMuon")
helper_0 = ROOT.SlotLocalOversampledHistogram1D(9, "h_oversampled", "Oversampled nMuon Histogram", 5, -0.5, 4.5)
h_oversampled = rdf.Book(helper_0, ("event", "nMuon"))
helper_1 = ROOT.SlotLocalOversampledHistogram1D(9, "h_gen_over", "Oversampled Gen Histogram", 5, -0.5, 4.5)
h_gen_over = rdf.Book(helper_1, ("event", "nGenJet"))
helper_2 = ROOT.SlotLocalOversampledHistogram1D(9, "h_sel_over", "Oversampled nMuon (nJet >= 30) >= 2", 5, -0.5, 4.5)
h_sel_over = rdf_sel.Book(helper_2, ("event", "nMuon"))
h_gen_base = rdf_base.Histo1D(("h_gen_base_1", "Base Gen nJet Histogram;Number of Gen Jets;Entries", 5, -0.5, 4.5), "nGenJet")
h_sel_base = rdf_base_sel.Histo1D(("h_sel_base", "Base nMuon (nJet >= 30) >= 2;Number of Muons;Entries", 5, -0.5, 4.5), "nMuon")
print("== nMuon ==")
h_oversampled.Print("all")
h_base.Sumw2()
h_base.Print("all")
print("\n")
print("== nGenJet ==")
h_gen_over.Print("all")
h_gen_base.Sumw2()
h_gen_base.Print("all")
print("\n== nMuon (nJet >= 30) >= 2 ==")
h_sel_over.Print("all")
h_sel_base.Sumw2()
h_sel_base.Print("all")
