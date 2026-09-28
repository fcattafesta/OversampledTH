#pragma once

#include <TH1.h>
#include <TROOT.h>

#include <ROOT/RDF/RActionImpl.hxx>
#include <ROOT/RVec.hxx>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

class TTreeReader;

namespace oversampled_detail {

  using EventId = unsigned long;

  inline void ValidateFactor(int factor) {
    if (factor < 1) {
      throw std::invalid_argument("Oversampling factor must be at least 1");
    }
  }

  template <typename TH>
  std::shared_ptr<TH> MakeHistogram(std::string_view name, std::string_view title, int bins, double min, double max) {
    static_assert(std::is_base_of_v<TH1, TH>, "TH must derive from TH1");
    auto hist = std::make_shared<TH>(std::string(name).c_str(), std::string(title).c_str(), bins, min, max);
    hist->SetDirectory(nullptr);
    hist->Sumw2();
    return hist;
  }

  template <typename TH>
  std::shared_ptr<TH> CloneEmptyShared(const TH &source) {
    auto hist = std::shared_ptr<TH>(static_cast<TH *>(source.Clone()));
    hist->SetDirectory(nullptr);
    hist->Reset();
    return hist;
  }

  template <typename TH>
  std::unique_ptr<TH> CloneEmpty(const TH &source) {
    auto hist = std::unique_ptr<TH>(static_cast<TH *>(source.Clone()));
    hist->SetDirectory(nullptr);
    hist->Reset();
    return hist;
  }

  // One Fill per generated event and bin makes Sumw2 track event-level contributions.
  template <typename TH>
  void FillEvent(TH &output, const TH &event, int factor) {
    for (int bin = 0; bin <= event.GetNbinsX() + 1; ++bin) {
      const double content = event.GetBinContent(bin);
      if (content != 0.0) {
        output.Fill(event.GetBinCenter(bin), content / factor);
      }
    }
  }

  // Preserve the intermediate bin precision of the supported TH1F/TH1D types.
  template <typename TH>
  using BinContent = std::conditional_t<std::is_base_of_v<TH1F, TH>, float, double>;

  template <typename Content>
  struct SparseEvent {
    struct Bin {
      int index;
      Content content;
    };
    std::vector<Bin> bins;  // Sorted, unique bin indices; zero totals are omitted.
    std::size_t pieces = 0;

    void Merge(SparseEvent &&other) {
      pieces += other.pieces;
      if (bins.empty()) {
        bins = std::move(other.bins);
        return;
      }
      if (other.bins.empty()) {
        return;
      }
      std::vector<Bin> combined;
      combined.reserve(bins.size() + other.bins.size());
      auto left = bins.begin();
      auto right = other.bins.begin();
      while (left != bins.end() && right != other.bins.end()) {
        if (left->index < right->index) {
          combined.push_back(*left++);
        } else if (right->index < left->index) {
          combined.push_back(*right++);
        } else {
          const Content total = static_cast<Content>(left->content + right->content);
          if (total != Content{}) {
            combined.push_back({left->index, total});
          }
          ++left;
          ++right;
        }
      }
      combined.insert(combined.end(), left, bins.end());
      combined.insert(combined.end(), right, other.bins.end());
      bins = std::move(combined);
    }
  };

  template <typename TH>
  SparseEvent<BinContent<TH>> ExtractEvent(const TH &hist) {
    SparseEvent<BinContent<TH>> event;
    event.pieces = 1;
    for (int bin = 0; bin <= hist.GetNbinsX() + 1; ++bin) {
      const auto content = static_cast<BinContent<TH>>(hist.GetBinContent(bin));
      if (content != BinContent<TH>{}) {
        event.bins.push_back({bin, content});
      }
    }
    return event;
  }

  template <typename TH>
  void FillEvent(TH &output, const SparseEvent<BinContent<TH>> &event, int factor) {
    for (const auto &bin : event.bins) {
      output.Fill(output.GetBinCenter(bin.index), static_cast<double>(bin.content) / factor);
    }
  }

  template <typename T>
  struct IsRVec : std::false_type {};

  template <typename T>
  struct IsRVec<ROOT::VecOps::RVec<T>> : std::true_type {};

  template <typename TH, typename Value>
  void FillValues(TH &hist, const Value &value, double weight) {
    if constexpr (IsRVec<std::decay_t<Value>>::value) {
      for (const auto &item : value) {
        hist.Fill(item, weight);
      }
    } else {
      hist.Fill(value, weight);
    }
  }

  template <typename TH>
  struct SlotState {
    std::unique_ptr<TH> event;
    std::unique_ptr<TH> output;
    EventId current = 0;
    bool hasCurrent = false;

    explicit SlotState(const TH &model) : event(CloneEmpty(model)), output(CloneEmpty(model)) {}

    void Flush(int factor) {
      if (hasCurrent) {
        FillEvent(*output, *event, factor);
        event->Reset();
        hasCurrent = false;
      }
    }
  };

  template <typename TH, typename State>
  std::vector<std::unique_ptr<State>> MakeSlots(const TH &model, unsigned int count) {
    std::vector<std::unique_ptr<State>> slots;
    slots.reserve(count);
    for (unsigned int slot = 0; slot < count; ++slot) {
      slots.push_back(std::make_unique<State>(model));
    }
    return slots;
  }

  inline unsigned int SlotCount() { return ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1U; }

}  // namespace oversampled_detail

template <typename TH>
class SequentialOversampledHistogram : public ROOT::Detail::RDF::RActionImpl<SequentialOversampledHistogram<TH>> {
public:
  using Result_t = TH;

  SequentialOversampledHistogram(
      int factor, std::string_view name, std::string_view title, int bins, double min, double max)
      : fFactor(factor),
        fResult(oversampled_detail::MakeHistogram<TH>(name, title, bins, min, max)),
        fEvent(oversampled_detail::CloneEmpty(*fResult)) {
    oversampled_detail::ValidateFactor(factor);
    if (ROOT::IsImplicitMTEnabled()) {
      throw std::runtime_error("SequentialOversampledHistogram requires implicit MT to be disabled");
    }
  }

  SequentialOversampledHistogram(const SequentialOversampledHistogram &other)
      : fFactor(other.fFactor),
        fResult(oversampled_detail::CloneEmptyShared(*other.fResult)),
        fEvent(oversampled_detail::CloneEmpty(*fResult)) {}

  SequentialOversampledHistogram &operator=(const SequentialOversampledHistogram &other) {
    if (this != &other) {
      *this = SequentialOversampledHistogram(other);
    }
    return *this;
  }

  SequentialOversampledHistogram(SequentialOversampledHistogram &&) = default;
  SequentialOversampledHistogram &operator=(SequentialOversampledHistogram &&) = default;

  std::shared_ptr<TH> GetResultPtr() const { return fResult; }
  void Initialize() {}
  void InitTask(TTreeReader *, unsigned int) {}

  template <typename Value>
  void Exec(unsigned int, oversampled_detail::EventId eventId, const Value &value, double weight = 1.0) {
    if (fHasCurrent && eventId != fCurrent) {
      Flush();
    }
    fCurrent = eventId;
    fHasCurrent = true;
    oversampled_detail::FillValues(*fEvent, value, weight);
  }

  void Flush() {
    if (fHasCurrent) {
      oversampled_detail::FillEvent(*fResult, *fEvent, fFactor);
      fEvent->Reset();
      fHasCurrent = false;
    }
  }

  void Finalize() { Flush(); }
  std::string GetActionName() { return "SequentialOversampledHistogram"; }

private:
  int fFactor;
  std::shared_ptr<TH> fResult;
  std::unique_ptr<TH> fEvent;
  oversampled_detail::EventId fCurrent = 0;
  bool fHasCurrent = false;
};

template <typename TH>
class SlotLocalOversampledHistogram : public ROOT::Detail::RDF::RActionImpl<SlotLocalOversampledHistogram<TH>> {
public:
  using Result_t = TH;

  SlotLocalOversampledHistogram(
      int factor, std::string_view name, std::string_view title, int bins, double min, double max)
      : fFactor(factor),
        fResult(oversampled_detail::MakeHistogram<TH>(name, title, bins, min, max)),
        fSlots(oversampled_detail::MakeSlots<TH, State>(*fResult, oversampled_detail::SlotCount())) {
    oversampled_detail::ValidateFactor(factor);
  }

  SlotLocalOversampledHistogram(const SlotLocalOversampledHistogram &other)
      : fFactor(other.fFactor),
        fResult(oversampled_detail::CloneEmptyShared(*other.fResult)),
        fSlots(oversampled_detail::MakeSlots<TH, State>(*fResult, other.fSlots.size())) {}

  SlotLocalOversampledHistogram &operator=(const SlotLocalOversampledHistogram &other) {
    if (this != &other) {
      *this = SlotLocalOversampledHistogram(other);
    }
    return *this;
  }

  SlotLocalOversampledHistogram(SlotLocalOversampledHistogram &&) = default;
  SlotLocalOversampledHistogram &operator=(SlotLocalOversampledHistogram &&) = default;

  std::shared_ptr<TH> GetResultPtr() const { return fResult; }
  void Initialize() {}
  void InitTask(TTreeReader *, unsigned int) {}

  template <typename Value>
  void Exec(unsigned int slot, oversampled_detail::EventId eventId, const Value &value, double weight = 1.0) {
    auto &state = *fSlots.at(slot);
    if (!state.hasCurrent || eventId != state.current) {
      if (!state.seen.insert(eventId).second) {
        state.revisited.insert(eventId);
      }
      if (state.hasCurrent) {
        state.Flush(fFactor);
      }
    }
    state.current = eventId;
    state.hasCurrent = true;
    oversampled_detail::FillValues(*state.event, value, weight);
  }

  void Flush() {
    for (auto &state : fSlots) {
      state->Flush(fFactor);
      fResult->Add(state->output.get());
      state->output->Reset();
    }
  }

  void Finalize() {
    Flush();
    std::unordered_map<oversampled_detail::EventId, unsigned int> counts;
    for (const auto &state : fSlots) {
      for (const auto id : state->seen) {
        ++counts[id];
      }
    }
    unsigned int split = 0;
    for (const auto &[id, count] : counts) {
      (void)id;
      split += count > 1;
    }
    std::unordered_set<oversampled_detail::EventId> revisited;
    for (const auto &state : fSlots) {
      revisited.insert(state->revisited.begin(), state->revisited.end());
    }
    std::cout << "SlotLocalOversampledHistogram diagnostics: unique genEvents=" << counts.size()
              << ", split across slots=" << split << ", revisited within slots=" << revisited.size() << '\n';
  }
  std::string GetActionName() { return "SlotLocalOversampledHistogram"; }

private:
  struct State : oversampled_detail::SlotState<TH> {
    using oversampled_detail::SlotState<TH>::SlotState;
    std::unordered_set<oversampled_detail::EventId> seen;
    std::unordered_set<oversampled_detail::EventId> revisited;
  };

  int fFactor;
  std::shared_ptr<TH> fResult;
  std::vector<std::unique_ptr<State>> fSlots;
};

template <typename TH>
class RangeAwareOversampledHistogram : public ROOT::Detail::RDF::RActionImpl<RangeAwareOversampledHistogram<TH>> {
public:
  using Result_t = TH;

  RangeAwareOversampledHistogram(
      int factor, std::string_view name, std::string_view title, int bins, double min, double max)
      : fFactor(factor),
        fResult(oversampled_detail::MakeHistogram<TH>(name, title, bins, min, max)),
        fSlots(oversampled_detail::MakeSlots<TH, State>(*fResult, oversampled_detail::SlotCount())) {
    oversampled_detail::ValidateFactor(factor);
  }

  RangeAwareOversampledHistogram(const RangeAwareOversampledHistogram &other)
      : fFactor(other.fFactor),
        fResult(oversampled_detail::CloneEmptyShared(*other.fResult)),
        fSlots(oversampled_detail::MakeSlots<TH, State>(*fResult, other.fSlots.size())) {}

  RangeAwareOversampledHistogram &operator=(const RangeAwareOversampledHistogram &other) {
    if (this != &other) {
      *this = RangeAwareOversampledHistogram(other);
    }
    return *this;
  }

  RangeAwareOversampledHistogram(RangeAwareOversampledHistogram &&) = default;
  RangeAwareOversampledHistogram &operator=(RangeAwareOversampledHistogram &&) = default;

  std::shared_ptr<TH> GetResultPtr() const { return fResult; }
  void Initialize() {}
  void InitTask(TTreeReader *, unsigned int) {}

  template <typename Value>
  void Exec(unsigned int slot,
            oversampled_detail::EventId eventId,
            const Value &value,
            unsigned long rangeBegin,
            unsigned long rangeEnd,
            double weight = 1.0) {
    auto &state = *fSlots.at(slot);
    if (!state.hasRange || state.rangeBegin != rangeBegin || state.rangeEnd != rangeEnd) {
      CloseRange(state);
      state.hasRange = true;
      state.rangeBegin = rangeBegin;
      state.rangeEnd = rangeEnd;
      state.hasFirst = false;
    }
    if (state.hasCurrent && eventId != state.current) {
      if (state.current == state.first) {
        DeferBoundary(state);
      } else {
        state.Flush(fFactor);
      }
    }
    if (!state.hasCurrent) {
      state.current = eventId;
      state.hasCurrent = true;
      if (!state.hasFirst) {
        state.first = eventId;
        state.hasFirst = true;
      }
    }
    oversampled_detail::FillValues(*state.event, value, weight);
  }

  void Flush() {
    for (auto &state : fSlots) {
      CloseRange(*state);
      fResult->Add(state->output.get());
      state->output->Reset();
    }

    BoundaryMap merged;
    for (auto &state : fSlots) {
      while (!state->boundary.empty()) {
        auto node = state->boundary.extract(state->boundary.begin());
        auto target = merged.find(node.key());
        if (target == merged.end()) {
          // Transfer the map node and its vector without cloning or copying bins.
          merged.insert(std::move(node));
        } else {
          target->second.Merge(std::move(node.mapped()));
        }
      }
      BoundaryMap{}.swap(state->boundary);  // Release consumed slot bucket storage too.
    }
    unsigned int split = 0;
    for (const auto &[id, event] : merged) {
      (void)id;
      oversampled_detail::FillEvent(*fResult, event, fFactor);
      split += event.pieces > 1;
    }
    std::cout << "RangeAwareOversampledHistogram diagnostics: boundary genEvents=" << merged.size()
              << ", split across ranges=" << split << '\n';
  }

  void Finalize() { Flush(); }
  std::string GetActionName() { return "RangeAwareOversampledHistogram"; }

private:
  using BoundaryEvent = oversampled_detail::SparseEvent<oversampled_detail::BinContent<TH>>;
  using BoundaryMap = std::unordered_map<oversampled_detail::EventId, BoundaryEvent>;

  struct State : oversampled_detail::SlotState<TH> {
    using oversampled_detail::SlotState<TH>::SlotState;
    bool hasRange = false;
    unsigned long rangeBegin = 0;
    unsigned long rangeEnd = 0;
    bool hasFirst = false;
    oversampled_detail::EventId first = 0;
    BoundaryMap boundary;
  };

  void DeferBoundary(State &state) {
    state.boundary[state.current].Merge(oversampled_detail::ExtractEvent(*state.event));
    state.event->Reset();
    state.hasCurrent = false;
  }

  void CloseRange(State &state) {
    if (state.hasRange && state.hasCurrent) {
      DeferBoundary(state);
    }
    state.hasRange = false;
    state.hasFirst = false;
  }

  int fFactor;
  std::shared_ptr<TH> fResult;
  std::vector<std::unique_ptr<State>> fSlots;
};

using SequentialOversampledHistogram1D = SequentialOversampledHistogram<TH1D>;
using SequentialOversampledHistogram1F = SequentialOversampledHistogram<TH1F>;
using SlotLocalOversampledHistogram1D = SlotLocalOversampledHistogram<TH1D>;
using SlotLocalOversampledHistogram1F = SlotLocalOversampledHistogram<TH1F>;
using RangeAwareOversampledHistogram1D = RangeAwareOversampledHistogram<TH1D>;
using RangeAwareOversampledHistogram1F = RangeAwareOversampledHistogram<TH1F>;
