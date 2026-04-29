import argparse
import os
import ROOT

include_path = os.path.join(os.path.dirname(__file__), os.pardir, "include")
ROOT.gInterpreter.AddIncludePath(include_path)
ROOT.gInterpreter.Declare('#include "DumbOversampledTH.h"')

parser = argparse.ArgumentParser(description="Test the DumbOversampledTH class.")
parser.add_argument("-j", type=int, default=1, help="Number of threads to use (default: 1)")
args = parser.parse_args()

if args.j > 1:
    print(f"Starting test with implicit multi-threading enabled. Using {args.j} threads.")
    ROOT.EnableImplicitMT(args.j)
else:
    print("Starting test with implicit multi-threading disabled.")

root_file = "/scratchnvme/store/mc/RunIII2024Summer24NanoAODv15/DY2Mu_2Jets_MLL-105to160-April2026_IreneFakes_Oversampling9_FlashSim/260421_142939/0000/tree_10.root"
file = ROOT.TFile.Open(root_file)
events = file.Get("Events")
print(f"Number of events: {events.GetEntries()}")

rdf = ROOT.RDataFrame("Events", root_file).Define("nGenJet", "static_cast<unsigned int>(GenJet_pt.size())")
rdf_base = rdf.Filter("fold == 0")
rdf_sel = rdf.Filter("ROOT::VecOps::Sum(Jet_pt > 30) >= 2")
rdf_base_sel = rdf_base.Filter("ROOT::VecOps::Sum(Jet_pt > 30) >= 2")

h_base = rdf_base.Histo1D(("h_base_nMuon", "Base nMuon Histogram;Number of Muons;Entries", 5, -0.5, 4.5), "nMuon")
helper_0 = ROOT.DumbOversampledTH1D(9, "h_oversampled", "Oversampled nMuon Histogram", 5, -0.5, 4.5)
h_oversampled = rdf.Book(helper_0, ("event", "nMuon"))
helper_1 = ROOT.DumbOversampledTH1D(9, "h_gen_over", "Oversampled Gen Histogram", 5, -0.5, 4.5)
h_gen_over = rdf.Book(helper_1, ("event", "nGenJet"))
helper_2 = ROOT.DumbOversampledTH1D(9, "h_sel_over", "Oversampled nMuon (nJet >= 30) >= 2", 5, -0.5, 4.5)
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