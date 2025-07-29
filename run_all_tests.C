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

  // Create jackknife object but keep a pointer to access methods later
  std::cout<< "----- JK book ----"<<std::endl;
  auto jackknifeTH = std::make_shared<OversampledTH<TH1F>>("OSHist_jackknife", "OversampledTH Jackknife", nBins, xMin, xMax, oversamplingFactor, 100);
  auto result_jackknife = rdf.Book(*jackknifeTH, {"event", histName.c_str()});
  std::cout<< "----- JK book done ----"<<std::endl;

  PrototypeOversampledTH<TH1F> firstProto{"FirstProto", "FirstPrototype", nBins, xMin, xMax};
  auto result_first = rdf.Book(std::move(firstProto), {"event", histName.c_str()});

  auto regular_result = hist_regular.GetValue();
  auto oversampled_result = result_oversampled.GetValue();
  auto first_result = result_first.GetValue();

  std::cout <<"--- JK GetValue ----"<<std::endl;
  auto jackknife_result = result_jackknife.GetValue();
  std::cout <<"--- end JK ----"<<std::endl;

  std::cout << "Regular Histogram - Entries: " << regular_result.GetEntries() << ", Mean: " << regular_result.GetMean()
            << std::endl;
  std::cout << "OversampledTH - Entries: " << oversampled_result.GetEntries()
            << ", Mean: " << oversampled_result.GetMean() << std::endl;
  std::cout << "FirstPrototype - Entries: " << first_result.GetEntries() << ", Mean: " << first_result.GetMean()
            << std::endl;
  std::cout << "Jackknife - Entries: " << jackknife_result.GetEntries() << ", Mean: " << jackknife_result.GetMean()
            << std::endl;

  auto outputFile = TFile::Open("comparison.root", "RECREATE");
  regular_result.Write("regular");
  oversampled_result.Write("oversampled");
  first_result.Write("first_prototype");
  jackknife_result.Write("jackknife");
  
  // Now we can access the covariance matrix through the shared pointer
  std::cout<<"---- JK COV -----" <<std::endl;
  auto covJackknife = jackknifeTH->getCovJackKnife();
  covJackknife->Write("cov_jackknife");
  auto avgJackknife = jackknifeTH->getJackknifeAverage();
  avgJackknife->Write("avg_jackkinfe");
  std::cout<<"---- JK COV done -----" <<std::endl;

  if (true){
    // Plot and save jackknife result and covariance matrix as PDF
    TCanvas c1("c1", "Jackknife Results", 1200, 600);
    c1.Divide(2,1);

    // Plot jackknife result histogram
    c1.cd(1);
    // set the errorbars from the diagonal element of the covariance matrix
    avgJackknife->SetStats(kFALSE);
    for(int i=1 ;i<=avgJackknife->GetNbinsX();++i )
    {
      avgJackknife->SetBinError(i, TMath::Sqrt( covJackknife->GetBinContent(i,i) ));
    }
    avgJackknife->SetLineColor(kBlue);
    avgJackknife->SetTitle("Jackknife Result");
    avgJackknife->Scale(1.0 / oversamplingFactor); // Scale by oversampling factor as others?
    avgJackknife->Draw("E1");
    // 

    oversampled_result.SetLineColor(kRed);
    oversampled_result.Draw("SAME E1");

    // Plot jackknife covariance matrix
    c1.cd(2);
    covJackknife->SetTitle("Jackknife Covariance Matrix");
    covJackknife->Draw("COLZ");

    // Save to PDF
    c1.SaveAs("jackknife_results.pdf");

    }

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
