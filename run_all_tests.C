#include <TFile.h>
#include <TH1F.h>

#include <ROOT/RDataFrame.hxx>
#include <chrono>
#include <iostream>

#include "include/FirstPrototype.h"
#include "include/OversampledTH.h"

double timeFunction(std::function<void()> func) {
  auto start = std::chrono::high_resolution_clock::now();
  func();
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  return duration.count();
}

int test_comparison() {
  std::cout << "=== Test: Histogram Comparison ===\n";

  auto rdf = ROOT::RDataFrame("Events", "tree_7.root")
                 .Define("Jet_p4",
                         "ROOT::VecOps::Construct<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>>(Jet_pt, "
                         "Jet_eta, Jet_phi, Jet_mass)")
                 .Filter("nJet >= 2")
                 .Define("Dijet_mass", "(Jet_p4[0] + Jet_p4[1]).M()");
  auto rdf_wo_os = rdf.Filter("fold == 0");

  std::string histName = "Dijet_mass";
  int nBins = 100;
  double xMin = 0.0;
  double xMax = 300.0;
  float oversamplingFactor = 10.0;

  auto hist_regular = rdf_wo_os.Histo1D({"regular", "Regular Histogram", nBins, xMin, xMax}, histName.c_str());

  OversampledTH<TH1F> oversampledTH{"OSHist", "OversampledTH", nBins, xMin, xMax, oversamplingFactor};
  auto result_oversampled = rdf.Book(std::move(oversampledTH), {"event", histName.c_str()});

  PrototypeOversampledTH<TH1F> firstProto{"FirstProto", "FirstPrototype", nBins, xMin, xMax};
  auto result_first = rdf.Book(std::move(firstProto), {"event", histName.c_str()});

  auto regular_result = hist_regular.GetValue();
  auto oversampled_result = result_oversampled.GetValue();
  auto first_result = result_first.GetValue();

  std::cout << "Regular Histogram - Entries: " << regular_result.GetEntries() << ", Mean: " << regular_result.GetMean()
            << std::endl;
  std::cout << "OversampledTH - Entries: " << oversampled_result.GetEntries()
            << ", Mean: " << oversampled_result.GetMean() << std::endl;
  std::cout << "FirstPrototype - Entries: " << first_result.GetEntries() << ", Mean: " << first_result.GetMean()
            << std::endl;

  auto outputFile = TFile::Open("comparison.root", "RECREATE");
  regular_result.Write("regular");
  oversampled_result.Write("oversampled");
  first_result.Write("first_prototype");
  outputFile->Close();

  return 0;
}

int test_sequential_vs_concurrent() {
  std::cout << "\n=== Test: Sequential vs Concurrent ===\n";

  std::string histName = "Dijet_mass";
  int nBins = 100;
  double xMin = 0.0;
  double xMax = 300.0;
  float oversamplingFactor = 10.0;

  // Sequential
  ROOT::DisableImplicitMT();
  std::cout << "Sequential mode:\n";
  double time_seq = timeFunction([&]() {
    auto rdf = ROOT::RDataFrame("Events", "tree_7.root");
    OversampledTH<TH1F> histHelper{"OSHist_seq", "Sequential", nBins, xMin, xMax, oversamplingFactor};
    auto result = rdf.Book(std::move(histHelper), {"event", histName.c_str()});
    auto final_result = result.GetValue();
  });
  std::cout << "  Time: " << time_seq << " ms\n";

  // Concurrent
  ROOT::EnableImplicitMT(5);
  std::cout << "Concurrent mode (5 threads):\n";
  double time_mt = timeFunction([&]() {
    auto rdf = ROOT::RDataFrame("Events", "tree_7.root");
    OversampledTH<TH1F> histHelper{"OSHist_mt", "Concurrent", nBins, xMin, xMax, oversamplingFactor};
    auto result = rdf.Book(std::move(histHelper), {"event", histName.c_str()});
    auto final_result = result.GetValue();
  });
  std::cout << "  Time: " << time_mt << " ms\n";
  std::cout << "  Speedup: " << time_seq / time_mt << "x\n";

  return 0;
}

int run_all_tests() {
  std::cout << "Running all tests...\n";

  test_comparison();
  // test_sequential_vs_concurrent();

  std::cout << "\nAll tests completed!\n";
  return 0;
}