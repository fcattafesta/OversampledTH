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

class TTreeReader;

template <typename TH>
class STOversampledTH : public ROOT::Detail::RDF::RActionImpl<STOversampledTH<TH>> {
public:
  using Result_t = TH;

private:
  int oversamplingFactor = 1;
  unsigned long fCurrentGenEvent = 0;
  bool fHasCurrentEvent = false;
  std::unique_ptr<TH> fCurrentEventHist;
  std::shared_ptr<TH> fFinalHist;

  std::unique_ptr<TH> MakeEmptyClone_() const {
    auto h = std::unique_ptr<TH>(static_cast<TH *>(fFinalHist->Clone()));
    h->SetDirectory(nullptr);
    h->Reset();
    return h;
  }

  void ResetState_() {
    fHasCurrentEvent = false;
    fCurrentGenEvent = 0;
    if (!fCurrentEventHist && fFinalHist) {
      fCurrentEventHist = MakeEmptyClone_();
    }
    if (fCurrentEventHist) {
      fCurrentEventHist->Reset();
    }
  }

  void FlushCurrentEvent_() {
    if (!fFinalHist || !fCurrentEventHist || !fHasCurrentEvent) {
      return;
    }
    const auto &hist = *fCurrentEventHist;
    for (int bin = 0; bin <= hist.GetNbinsX() + 1; ++bin) {
      const double content = hist.GetBinContent(bin);
      if (content == 0.0) {
        continue;
      }
      fFinalHist->Fill(hist.GetBinCenter(bin), content / static_cast<double>(oversamplingFactor));
    }
    fCurrentEventHist->Reset();
    fHasCurrentEvent = false;
  }

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
    fCurrentEventHist = MakeEmptyClone_();
  }
  STOversampledTH(const STOversampledTH &other)
      : oversamplingFactor(other.oversamplingFactor),
        fFinalHist(other.fFinalHist ? std::shared_ptr<TH>(static_cast<TH *>(other.fFinalHist->Clone())) : nullptr) {
    if (fFinalHist) {
      fFinalHist->SetDirectory(nullptr);
      fFinalHist->Reset();
      fCurrentEventHist = MakeEmptyClone_();
    }
    ResetState_();
  }
  STOversampledTH &operator=(const STOversampledTH &other) {
    if (this == &other) {
      return *this;
    }

    oversamplingFactor = other.oversamplingFactor;
    fFinalHist = other.fFinalHist ? std::shared_ptr<TH>(static_cast<TH *>(other.fFinalHist->Clone())) : nullptr;
    if (fFinalHist) {
      fFinalHist->SetDirectory(nullptr);
      fFinalHist->Reset();
      fCurrentEventHist = MakeEmptyClone_();
    } else {
      fCurrentEventHist.reset();
    }
    ResetState_();
    return *this;
  }
  STOversampledTH(STOversampledTH &&) = default;
  STOversampledTH &operator=(STOversampledTH &&) = default;
  ~STOversampledTH() = default;
  std::shared_ptr<TH> GetResultPtr() const { return fFinalHist; }
  void Initialize() {}
  void InitTask(TTreeReader *, unsigned int) {
    if (!fCurrentEventHist && fFinalHist) {
      fCurrentEventHist = MakeEmptyClone_();
    }
  }

  //   template <typename T>
  //   void Exec(unsigned long genEvent, ROOT::VecOps::RVec<T> values, double weight = 1);

  template <typename T>
  void Exec(unsigned int slot, unsigned long genEvent, T values, double weight = 1) {
    (void)slot;
    if (!fFinalHist) {
      std::cerr << "Error: fFinalHist is null" << std::endl;
      return;
    }
    if (!fCurrentEventHist) {
      fCurrentEventHist = MakeEmptyClone_();
    }

    if (!fHasCurrentEvent) {
      fHasCurrentEvent = true;
      fCurrentGenEvent = genEvent;
    } else if (genEvent != fCurrentGenEvent) {
      FlushCurrentEvent_();
      fHasCurrentEvent = true;
      fCurrentGenEvent = genEvent;
    }

    fCurrentEventHist->Fill(values, weight);
  }

  void Flush() { FlushCurrentEvent_(); }

  void Finalize() {
    //   std::cout << "Finalizing STOversampledTH." << std::endl;
    Flush();
  }
  std::string GetActionName() { return "STOversampledTH"; }
};

using STOversampledTH1D = STOversampledTH<TH1D>;
using STOversampledTH1F = STOversampledTH<TH1F>;