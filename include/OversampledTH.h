#pragma once

#include <TH1.h>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

template <typename TH>
class OversampledTH : public ROOT::Detail::RDF::RActionImpl<OversampledTH<TH>> {
public:
  using Result_t = TH;

private:
  unsigned int nSlots_;
  const float oversamplingFactor_{1.0f};
  long int lastFlush_ = -1;
  // { eventId -> (slot -> histogram) }
  std::unordered_map<long int, std::unordered_map<unsigned int, TH *>> partialHists_;
  // {slot -> current event id}
  std::vector<long int> currentEvents_;
  std::shared_ptr<TH> finalHistogram_;

public:
  OversampledTH(std::string_view name,
                std::string_view title,
                int nbin,
                double xmin,
                double xmax,
                float oversamplingFactor = 1.0f)
      : nSlots_(ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1),
        currentEvents_(nSlots_, -1),
        oversamplingFactor_(oversamplingFactor),
        finalHistogram_(std::make_shared<TH>(std::string(name).c_str(), std::string(title).c_str(), nbin, xmin, xmax)) {
  }
  // TODO: Implement copy, move constructors and assignment operators

  // Mandatory Interface

  std::shared_ptr<TH> GetResultPtr() const { return finalHistogram_; }

  void Initialize() {}

  void InitTask(TTreeReader *, unsigned int) {}

  template <typename ValueType>
  void Exec(unsigned int slot, long int event, ValueType values, double weight = 1) {
    if constexpr (is_rvec_v<ValueType>) {
      for (const auto &v : values) {
        exec_(slot, event, v, weight);
      }
    } else {
      exec_(slot, event, values, weight);
    }
  }
  void Finalize() { flush_(true); }
  std::string GetActionName() { return "OversampledTH"; }

  // TODO: Optional Interface

private:
  // Helper to check if a type is ROOT::VecOps::RVec<T>
  template <typename T>
  struct is_rvec : std::false_type {};
  template <typename T>
  struct is_rvec<ROOT::VecOps::RVec<T>> : std::true_type {};
  template <typename T>
  static constexpr bool is_rvec_v = is_rvec<T>::value;

  void exec_(unsigned int slot, long int event, double value, double weight = 1) {
    // Check if the event already has a histogram for this slot
    if (partialHists_[event].find(slot) == partialHists_[event].end()) {
      TH *h = (TH *)finalHistogram_->Clone();
      partialHists_[event].emplace(slot, h);
      partialHists_[event][slot]->Reset();
    }
    // Fill the histogram for the current event and slot
    partialHists_[event][slot]->Fill(value, weight);
    // If the event is different from the last recorded event for this slot, update and flush
    if (event != currentEvents_[slot]) {
      currentEvents_[slot] = event;
      flush_();
    }
  }

  void fillPartialHists_(const std::unordered_map<unsigned int, TH *> &partialHists) {
    // Fill the final histogram with the contents of the partial histograms
    for (const auto &[slot, histogram] : partialHists) {
      for (int bin = 0; bin <= histogram->GetNbinsX(); bin++) {
        finalHistogram_->Fill(histogram->GetBinCenter(bin), histogram->GetBinContent(bin) / oversamplingFactor_);
      }
      delete (histogram);
    }
  }

  void flush_(bool all = false) {
    // Check if we need to flush based on the minimum event ID
    auto minimumEventId = *std::min_element(currentEvents_.begin(), currentEvents_.end());
    // If the last flush was before the minimum event ID or if 'all' is true, flush the partial histograms
    if ((lastFlush_ < (minimumEventId - 1)) || all) {
      for (auto it = partialHists_.begin(); it != partialHists_.end();) {
        auto eventId = it->first;
        // Flush all the histograms for events with IDs less than the minimum event ID
        if (eventId < minimumEventId) {
          fillPartialHists_(it->second);
          lastFlush_ = eventId;
          it = partialHists_.erase(it);
        } else {
          ++it;
        }
      }
    }
  }
};
