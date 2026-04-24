/*
Multi-threaded oversampled histogram action for ROOT RDataFrame. 
Dumb because we don't care about same genEvent in different threads.
*/
#pragma once

#include <TH1.h>
#include <TROOT.h>

#include <ROOT/RDF/RActionImpl.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class TTreeReader;

template <typename TH>
class DumbOversampledTH : public ROOT::Detail::RDF::RActionImpl<DumbOversampledTH<TH>> {
public:
  using Result_t = TH;

private:
  int oversamplingFactor = 1;
  unsigned int nSlots = 1;
  std::vector<unsigned long> fCurrentGenEvent;
  std::vector<bool> fHasCurrentEvent;
  std::vector<std::unique_ptr<TH>> fCurrentEventHists;
  std::vector<std::unique_ptr<TH>> fSlotHists;
  std::vector<std::unordered_set<unsigned long>> fSeenGenEventsPerSlot;
  std::shared_ptr<TH> fFinalHist;

  std::unique_ptr<TH> MakeEmptyClone_() const {
    auto h = std::unique_ptr<TH>(static_cast<TH *>(fFinalHist->Clone()));
    h->SetDirectory(nullptr);
    h->Reset();
    return h;
  }

  void ClearHists_() {
    for (auto &h : fSlotHists) {
      h.reset();
    }
  }

public:
  DumbOversampledTH(
      int oversamplingFactor, std::string_view name, std::string_view title, int nbin, double xmin, double xmax)
      : oversamplingFactor(oversamplingFactor),
        nSlots(ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1),
        fCurrentGenEvent(nSlots, std::numeric_limits<unsigned long>::max()),
        fHasCurrentEvent(nSlots, false),
        fCurrentEventHists(nSlots),
        fSlotHists(nSlots),
        fSeenGenEventsPerSlot(nSlots),
        fFinalHist(std::make_shared<TH>(std::string(name).c_str(), std::string(title).c_str(), nbin, xmin, xmax)) {
    if (oversamplingFactor < 1) {
      throw std::invalid_argument("Oversampling factor must be at least 1");
    }
    fFinalHist->SetDirectory(nullptr);
  }
  DumbOversampledTH(const DumbOversampledTH &other)
      : oversamplingFactor(other.oversamplingFactor),
        nSlots(other.nSlots),
        fCurrentGenEvent(other.nSlots, std::numeric_limits<unsigned long>::max()),
        fHasCurrentEvent(other.nSlots, false),
        fCurrentEventHists(other.nSlots),
        fSlotHists(other.nSlots),
        fSeenGenEventsPerSlot(other.nSlots),
        fFinalHist(other.fFinalHist ? std::shared_ptr<TH>(static_cast<TH *>(other.fFinalHist->Clone())) : nullptr) {
    if (fFinalHist) {
      fFinalHist->SetDirectory(nullptr);
      fFinalHist->Reset();
    }
  }
  DumbOversampledTH &operator=(const DumbOversampledTH &other) {
    if (this == &other) {
      return *this;
    }

    ClearHists_();
    oversamplingFactor = other.oversamplingFactor;
    nSlots = other.nSlots;
    fCurrentGenEvent.assign(nSlots, std::numeric_limits<unsigned long>::max());
    fHasCurrentEvent.assign(nSlots, false);
    fCurrentEventHists.assign(nSlots, nullptr);
    fSlotHists.assign(nSlots, nullptr);
    fSeenGenEventsPerSlot.assign(nSlots, {});
    fFinalHist = other.fFinalHist ? std::shared_ptr<TH>(static_cast<TH *>(other.fFinalHist->Clone())) : nullptr;
    if (fFinalHist) {
      fFinalHist->SetDirectory(nullptr);
      fFinalHist->Reset();
    }
    return *this;
  }
  DumbOversampledTH(DumbOversampledTH &&) = default;
  DumbOversampledTH &operator=(DumbOversampledTH &&) = default;
  ~DumbOversampledTH() { ClearHists_(); }
  std::shared_ptr<TH> GetResultPtr() const { return fFinalHist; }
  void Initialize() {}
  void InitTask(TTreeReader *, unsigned int slot) {
    if (slot >= fSlotHists.size()) {
      return;
    }
    if (!fCurrentEventHists[slot]) {
      fCurrentEventHists[slot] = MakeEmptyClone_();
    }
    if (!fSlotHists[slot]) {
      fSlotHists[slot] = MakeEmptyClone_();
    }
  }

  void FlushCurrentEvent_(unsigned int slot) {
    if (slot >= fCurrentEventHists.size() || !fCurrentEventHists[slot] || !fSlotHists[slot]) {
      return;
    }
    auto &eventHist = *fCurrentEventHists[slot];
    auto &slotHist = *fSlotHists[slot];
    for (int bin = 0; bin <= eventHist.GetNbinsX() + 1; ++bin) {
      const double content = eventHist.GetBinContent(bin);
      // if (content == 0.0) {
      //   continue;
      // }
      slotHist.Fill(eventHist.GetBinCenter(bin), content / static_cast<double>(oversamplingFactor));
    }
    eventHist.Reset();
  }

  template <typename T>
  void FillOne_(unsigned int slot, unsigned long genEvent, const T &value, double weight) {
    if (!fFinalHist || slot >= fSlotHists.size()) {
      return;
    }
    if (!fCurrentEventHists[slot]) {
      fCurrentEventHists[slot] = MakeEmptyClone_();
    }
    if (!fSlotHists[slot]) {
      fSlotHists[slot] = MakeEmptyClone_();
    }
    fSeenGenEventsPerSlot[slot].insert(genEvent);
    if (!fHasCurrentEvent[slot]) {
      fHasCurrentEvent[slot] = true;
      fCurrentGenEvent[slot] = genEvent;
    } else if (fCurrentGenEvent[slot] != genEvent) {
      FlushCurrentEvent_(slot);
      fCurrentGenEvent[slot] = genEvent;
    }
    fCurrentEventHists[slot]->Fill(value, weight);
  }

  template <typename T>
  struct is_rvec : std::false_type {};
  template <typename T>
  struct is_rvec<ROOT::VecOps::RVec<T>> : std::true_type {};

  template <typename T>
  void Exec(unsigned int slot, unsigned long genEvent, T values, double weight = 1) {
    if constexpr (is_rvec<T>::value) {
      for (const auto &v : values) {
        FillOne_(slot, genEvent, v, weight);
      }
    } else {
      FillOne_(slot, genEvent, values, weight);
    }
  }

  void Flush() {
    if (!fFinalHist) {
      return;
    }
    for (unsigned int slot = 0; slot < nSlots; ++slot) {
      if (fHasCurrentEvent[slot]) {
        FlushCurrentEvent_(slot);
      }
      fHasCurrentEvent[slot] = false;
      fCurrentGenEvent[slot] = std::numeric_limits<unsigned long>::max();
    }
    for (auto &h : fSlotHists) {
      if (h) {
        fFinalHist->Add(h.get());
        h->Reset();
      }
    }
  }

  void Finalize() {
    Flush();

    std::unordered_map<unsigned long, unsigned int> eventSlotCounts;
    for (const auto &seenInSlot : fSeenGenEventsPerSlot) {
      for (const auto genEvent : seenInSlot) {
        eventSlotCounts[genEvent] += 1;
      }
    }

    unsigned long splitAcrossSlots = 0;
    for (const auto &kv : eventSlotCounts) {
      if (kv.second > 1) {
        ++splitAcrossSlots;
      }
    }

    std::cout << "DumbOversampledTH diagnostics: unique genEvents=" << eventSlotCounts.size()
              << ", split across slots=" << splitAcrossSlots << std::endl;
  }
  std::string GetActionName() { return "DumbOversampledTH"; }
};

using DumbOversampledTH1D = DumbOversampledTH<TH1D>;
using DumbOversampledTH1F = DumbOversampledTH<TH1F>;