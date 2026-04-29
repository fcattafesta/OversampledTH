/*
Range-aware oversampled histogram action for ROOT RDataFrame.
Assumes identical genEvent values are contiguous in entry order.
Only range-boundary events are deferred; interior events are flushed early.
*/
#pragma once

#include <TH1.h>
#include <TROOT.h>

#include <ROOT/RDF/RActionImpl.hxx>
#include <ROOT/RVec.hxx>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

class TTreeReader;

template <typename TH>
class RangeOversampledTH : public ROOT::Detail::RDF::RActionImpl<RangeOversampledTH<TH>> {
public:
  using Result_t = TH;

private:
  int fOversamplingFactor = 1;
  unsigned int fNSlots = 1;

  std::vector<bool> fHasCurrentRange;
  std::vector<unsigned long> fRangeBegin;
  std::vector<unsigned long> fRangeEnd;

  std::vector<bool> fHasCurrentEvent;
  std::vector<unsigned long> fCurrentEvent;
  std::vector<bool> fHasFirstEventInRange;
  std::vector<unsigned long> fFirstEventInRange;

  std::vector<std::unique_ptr<TH>> fCurrentEventHists;
  std::vector<std::unique_ptr<TH>> fSlotHists;
  std::vector<std::unordered_map<unsigned long, std::unique_ptr<TH>>> fBoundaryHistsPerSlot;

  std::shared_ptr<TH> fFinalHist;

  std::unique_ptr<TH> MakeEmptyClone_() const {
    auto h = std::unique_ptr<TH>(static_cast<TH *>(fFinalHist->Clone()));
    h->SetDirectory(nullptr);
    h->Reset();
    return h;
  }

  void EnsureSlotBuffers_(unsigned int slot) {
    if (slot >= fNSlots) {
      return;
    }
    if (!fCurrentEventHists[slot]) {
      fCurrentEventHists[slot] = MakeEmptyClone_();
    }
    if (!fSlotHists[slot]) {
      fSlotHists[slot] = MakeEmptyClone_();
    }
  }

  void FillScaledInto_(TH &target, const TH &eventHist) {
    for (int bin = 0; bin <= eventHist.GetNbinsX() + 1; ++bin) {
      const double content = eventHist.GetBinContent(bin);
      if (content == 0.0) {
        continue;
      }
      target.Fill(eventHist.GetBinCenter(bin), content / static_cast<double>(fOversamplingFactor));
    }
  }

  void FlushCurrentEventToSlot_(unsigned int slot) {
    if (slot >= fNSlots || !fCurrentEventHists[slot] || !fSlotHists[slot]) {
      return;
    }
    FillScaledInto_(*fSlotHists[slot], *fCurrentEventHists[slot]);
    fCurrentEventHists[slot]->Reset();
  }

  void AddCurrentEventToBoundary_(unsigned int slot) {
    if (slot >= fNSlots || !fCurrentEventHists[slot]) {
      return;
    }

    const auto eventId = fCurrentEvent[slot];
    auto &boundaryMap = fBoundaryHistsPerSlot[slot];
    auto it = boundaryMap.find(eventId);
    if (it == boundaryMap.end()) {
      auto h = MakeEmptyClone_();
      h->Add(fCurrentEventHists[slot].get());
      boundaryMap.emplace(eventId, std::move(h));
    } else {
      it->second->Add(fCurrentEventHists[slot].get());
    }
    fCurrentEventHists[slot]->Reset();
  }

  void StartRange_(unsigned int slot, unsigned long rangeBegin, unsigned long rangeEnd) {
    fHasCurrentRange[slot] = true;
    fRangeBegin[slot] = rangeBegin;
    fRangeEnd[slot] = rangeEnd;
    fHasFirstEventInRange[slot] = false;
  }

  void CloseRange_(unsigned int slot) {
    if (slot >= fNSlots || !fHasCurrentRange[slot]) {
      return;
    }

    if (fHasCurrentEvent[slot]) {
      // Last event in the range is always a boundary candidate.
      AddCurrentEventToBoundary_(slot);
      fHasCurrentEvent[slot] = false;
    }

    fHasCurrentRange[slot] = false;
    fHasFirstEventInRange[slot] = false;
    fRangeBegin[slot] = 0;
    fRangeEnd[slot] = 0;
  }

  template <typename T>
  void FillOne_(unsigned int slot,
                unsigned long genEvent,
                const T &value,
                unsigned long rangeBegin,
                unsigned long rangeEnd,
                double weight) {
    if (!fFinalHist || slot >= fNSlots) {
      return;
    }

    EnsureSlotBuffers_(slot);

    if (!fHasCurrentRange[slot]) {
      StartRange_(slot, rangeBegin, rangeEnd);
    } else if (fRangeBegin[slot] != rangeBegin || fRangeEnd[slot] != rangeEnd) {
      // New processing range for this slot: close previous range first.
      CloseRange_(slot);
      StartRange_(slot, rangeBegin, rangeEnd);
    }

    if (!fHasCurrentEvent[slot]) {
      fHasCurrentEvent[slot] = true;
      fCurrentEvent[slot] = genEvent;
      if (!fHasFirstEventInRange[slot]) {
        fHasFirstEventInRange[slot] = true;
        fFirstEventInRange[slot] = genEvent;
      }
    } else if (fCurrentEvent[slot] != genEvent) {
      // On event transition within the same range:
      // - first event in range -> boundary buffer
      // - interior events -> early flush into slot histogram
      if (fHasFirstEventInRange[slot] && fCurrentEvent[slot] == fFirstEventInRange[slot]) {
        AddCurrentEventToBoundary_(slot);
      } else {
        FlushCurrentEventToSlot_(slot);
      }
      fCurrentEvent[slot] = genEvent;
    }

    fCurrentEventHists[slot]->Fill(value, weight);
  }

public:
  RangeOversampledTH(
      int oversamplingFactor, std::string_view name, std::string_view title, int nbin, double xmin, double xmax)
      : fOversamplingFactor(oversamplingFactor),
        fNSlots(ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1),
        fHasCurrentRange(fNSlots, false),
        fRangeBegin(fNSlots, 0),
        fRangeEnd(fNSlots, 0),
        fHasCurrentEvent(fNSlots, false),
        fCurrentEvent(fNSlots, std::numeric_limits<unsigned long>::max()),
        fHasFirstEventInRange(fNSlots, false),
        fFirstEventInRange(fNSlots, std::numeric_limits<unsigned long>::max()),
        fCurrentEventHists(fNSlots),
        fSlotHists(fNSlots),
        fBoundaryHistsPerSlot(fNSlots),
        fFinalHist(std::make_shared<TH>(std::string(name).c_str(), std::string(title).c_str(), nbin, xmin, xmax)) {
    if (fOversamplingFactor < 1) {
      throw std::invalid_argument("Oversampling factor must be at least 1");
    }
    fFinalHist->SetDirectory(nullptr);
  }

  RangeOversampledTH(const RangeOversampledTH &other)
      : fOversamplingFactor(other.fOversamplingFactor),
        fNSlots(other.fNSlots),
        fHasCurrentRange(other.fNSlots, false),
        fRangeBegin(other.fNSlots, 0),
        fRangeEnd(other.fNSlots, 0),
        fHasCurrentEvent(other.fNSlots, false),
        fCurrentEvent(other.fNSlots, std::numeric_limits<unsigned long>::max()),
        fHasFirstEventInRange(other.fNSlots, false),
        fFirstEventInRange(other.fNSlots, std::numeric_limits<unsigned long>::max()),
        fCurrentEventHists(other.fNSlots),
        fSlotHists(other.fNSlots),
        fBoundaryHistsPerSlot(other.fNSlots),
        fFinalHist(other.fFinalHist ? std::shared_ptr<TH>(static_cast<TH *>(other.fFinalHist->Clone())) : nullptr) {
    if (fFinalHist) {
      fFinalHist->SetDirectory(nullptr);
      fFinalHist->Reset();
    }
  }

  RangeOversampledTH &operator=(const RangeOversampledTH &other) {
    if (this == &other) {
      return *this;
    }

    fOversamplingFactor = other.fOversamplingFactor;
    fNSlots = other.fNSlots;

    fHasCurrentRange.assign(fNSlots, false);
    fRangeBegin.assign(fNSlots, 0);
    fRangeEnd.assign(fNSlots, 0);

    fHasCurrentEvent.assign(fNSlots, false);
    fCurrentEvent.assign(fNSlots, std::numeric_limits<unsigned long>::max());
    fHasFirstEventInRange.assign(fNSlots, false);
    fFirstEventInRange.assign(fNSlots, std::numeric_limits<unsigned long>::max());

    fCurrentEventHists.assign(fNSlots, nullptr);
    fSlotHists.assign(fNSlots, nullptr);
    fBoundaryHistsPerSlot.assign(fNSlots, {});

    fFinalHist = other.fFinalHist ? std::shared_ptr<TH>(static_cast<TH *>(other.fFinalHist->Clone())) : nullptr;
    if (fFinalHist) {
      fFinalHist->SetDirectory(nullptr);
      fFinalHist->Reset();
    }

    return *this;
  }

  RangeOversampledTH(RangeOversampledTH &&) = default;
  RangeOversampledTH &operator=(RangeOversampledTH &&) = default;
  ~RangeOversampledTH() = default;

  std::shared_ptr<TH> GetResultPtr() const { return fFinalHist; }
  void Initialize() {}

  void InitTask(TTreeReader *, unsigned int slot) {
    if (slot >= fNSlots) {
      return;
    }
    EnsureSlotBuffers_(slot);
  }

  template <typename T>
  void Exec(unsigned int slot,
            unsigned long genEvent,
            const ROOT::VecOps::RVec<T> &values,
            unsigned long rangeBegin,
            unsigned long rangeEnd,
            double weight = 1) {
    for (const auto &v : values) {
      FillOne_(slot, genEvent, v, rangeBegin, rangeEnd, weight);
    }
  }

  template <typename T>
  void Exec(unsigned int slot,
            unsigned long genEvent,
            T value,
            unsigned long rangeBegin,
            unsigned long rangeEnd,
            double weight = 1) {
    FillOne_(slot, genEvent, value, rangeBegin, rangeEnd, weight);
  }

  void Flush() {
    if (!fFinalHist) {
      return;
    }

    for (unsigned int slot = 0; slot < fNSlots; ++slot) {
      if (fHasCurrentRange[slot]) {
        CloseRange_(slot);
      }
    }

    // Early-flushed interior events are already in slot histograms.
    for (auto &h : fSlotHists) {
      if (h) {
        fFinalHist->Add(h.get());
        h->Reset();
      }
    }

    // Reconcile deferred boundary events across slots/ranges.
    std::unordered_map<unsigned long, std::unique_ptr<TH>> mergedBoundary;
    std::unordered_map<unsigned long, unsigned int> boundaryPieces;

    for (unsigned int slot = 0; slot < fNSlots; ++slot) {
      auto &slotBoundary = fBoundaryHistsPerSlot[slot];
      for (auto &kv : slotBoundary) {
        const auto eventId = kv.first;
        boundaryPieces[eventId] += 1;

        auto it = mergedBoundary.find(eventId);
        if (it == mergedBoundary.end()) {
          auto h = MakeEmptyClone_();
          h->Add(kv.second.get());
          mergedBoundary.emplace(eventId, std::move(h));
        } else {
          it->second->Add(kv.second.get());
        }
      }
      slotBoundary.clear();
    }

    for (const auto &kv : mergedBoundary) {
      FillScaledInto_(*fFinalHist, *kv.second);
    }

    unsigned long splitAcrossRanges = 0;
    for (const auto &kv : boundaryPieces) {
      if (kv.second > 1) {
        ++splitAcrossRanges;
      }
    }

    std::cout << "RangeOversampledTH diagnostics: boundary genEvents=" << mergedBoundary.size()
              << ", split across ranges=" << splitAcrossRanges << std::endl;
  }

  void Finalize() { Flush(); }
  std::string GetActionName() { return "RangeOversampledTH"; }
};

using RangeOversampledTH1D = RangeOversampledTH<TH1D>;
using RangeOversampledTH1F = RangeOversampledTH<TH1F>;
