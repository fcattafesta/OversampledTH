/*
Single-threaded oversampled histogram action for ROOT RDataFrame. 
*/
#pragma once

#include <TH1.h>
#include <TROOT.h>

#include <ROOT/RDF/RActionImpl.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>

class TTreeReader;

template <typename TH>
class STOversampledTH : public ROOT::Detail::RDF::RActionImpl<STOversampledTH<TH>> {
public:
  using Result_t = TH;

private:
  int oversamplingFactor = 1;
  long int lastGenEvent = -1;
  std::unordered_map<unsigned long, std::unique_ptr<TH>> fHists;
  std::shared_ptr<TH> fFinalHist;

  void ClearHists_() { fHists.clear(); }

public:
  STOversampledTH(
      int oversamplingFactor, std::string_view name, std::string_view title, int nbin, double xmin, double xmax)
      : oversamplingFactor(oversamplingFactor),
        fFinalHist(std::make_shared<TH>(std::string(name).c_str(), std::string(title).c_str(), nbin, xmin, xmax)) {
    if (oversamplingFactor < 1) {
      throw std::invalid_argument("Oversampling factor must be at least 1");
    }
    if (ROOT::IsImplicitMTEnabled()) {
      throw std::runtime_error(
          "STOversampledTH is designed for single-threaded execution. Please disable implicit multi-threading.");
    }
    fFinalHist->SetDirectory(nullptr);
  }
  STOversampledTH(const STOversampledTH &other)
      : oversamplingFactor(other.oversamplingFactor),
        lastGenEvent(-1),
        fFinalHist(other.fFinalHist ? std::shared_ptr<TH>(static_cast<TH *>(other.fFinalHist->Clone())) : nullptr) {
    if (fFinalHist) {
      fFinalHist->SetDirectory(nullptr);
    }
  }
  STOversampledTH &operator=(const STOversampledTH &other) {
    if (this == &other) {
      return *this;
    }

    ClearHists_();
    oversamplingFactor = other.oversamplingFactor;
    lastGenEvent = -1;
    fFinalHist = other.fFinalHist ? std::shared_ptr<TH>(static_cast<TH *>(other.fFinalHist->Clone())) : nullptr;
    if (fFinalHist) {
      fFinalHist->SetDirectory(nullptr);
    }
    return *this;
  }
  STOversampledTH(STOversampledTH &&) = default;
  STOversampledTH &operator=(STOversampledTH &&) = default;
  ~STOversampledTH() { ClearHists_(); }
  std::shared_ptr<TH> GetResultPtr() const { return fFinalHist; }
  void Initialize() {}
  void InitTask(TTreeReader *, unsigned int) {}

  //   template <typename T>
  //   void Exec(unsigned long genEvent, ROOT::VecOps::RVec<T> values, double weight = 1);

  template <typename T>
  void Exec(unsigned int slot, unsigned long genEvent, T values, double weight = 1) {
    if (!fFinalHist) {
      std::cerr << "Error: fFinalHist is null" << std::endl;
      return;
    }

    if (genEvent != lastGenEvent) {
    //   std::cout << "GenEvent changed from " << lastGenEvent << " to " << genEvent << std::endl;
    //   std::cout << "Flushing histograms @ genEvent: " << lastGenEvent << std::endl;
      Flush();
      lastGenEvent = genEvent;
    }

    if (fHists.find(genEvent) == fHists.end()) {
      auto h = std::unique_ptr<TH>(static_cast<TH *>(fFinalHist->Clone()));
      h->SetDirectory(nullptr);
      h->Reset();
      fHists.emplace(genEvent, std::move(h));
    }
    fHists[genEvent]->Fill(values, weight);
  }

  void Flush() {
    for (const auto &kv : fHists) {
      const auto &hist = *kv.second;
      const auto &genEvent = kv.first;
      //   std::cout << "     -> Flushing histogram for genEvent: " << genEvent << std::endl;
      for (size_t bin = 0; bin <= hist.GetNbinsX() + 1; bin++) {
        fFinalHist->Fill(hist.GetBinCenter(bin), hist.GetBinContent(bin) / static_cast<double>(oversamplingFactor));
      }
    }
    fHists.clear();
  }

  void Finalize() {
    //   std::cout << "Finalizing STOversampledTH." << std::endl;
    Flush();
  }
  std::string GetActionName() { return "STOversampledTH"; }
};

using STOversampledTH1D = STOversampledTH<TH1D>;
using STOversampledTH1F = STOversampledTH<TH1F>;