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
  using Result2D_t = TH2D; // Assuming Result2_t is a 2D histogram, adjust as needed

private:
  unsigned int nSlots_;
  const float oversamplingFactor_{1.0f};
  long int lastFlush_ = -1;
  // { eventId -> (slot -> histogram) }
  std::unordered_map<long int, std::unordered_map<unsigned int, TH *>> partialHists_;
  // {slot -> current event id}
  std::vector<long int> currentEvents_;
  std::shared_ptr<TH> finalHistogram_;

  //introduce a vector of Histograms for jackknife resampling
  int J = 0; // number of jacknife resampling
  std::vector<std::shared_ptr<TH>> jackknifeHistsVec_;
  std::vector<long> jackknifeN_; // number of event in each jackknife resampling hist
  bool hasCovJackKnife_={false};
  std::shared_ptr<Result2D_t> jackknifeCov_;
  std::shared_ptr<TH> jackknifeAverage_;

public:
  OversampledTH(std::string_view name,
                std::string_view title,
                int nbin,
                double xmin,
                double xmax,
                float oversamplingFactor = 1.0f,
                int jackknifeResampling = 0 // if 0 no jackknife resampling, otherwise the number of d-deleted resampling 
              )
      : nSlots_(ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1),
        currentEvents_(nSlots_, -1),
        oversamplingFactor_(oversamplingFactor),
        J (jackknifeResampling),
        finalHistogram_(std::make_shared<TH>(std::string(name).c_str(), std::string(title).c_str(), nbin, xmin, xmax)) {

          if (J>1){ // init jack knife resampling vector
            jackknifeHistsVec_.reserve(J);
            jackknifeN_.resize(J, 0); // initialize the number of events in each jackknife resampling
            for (int i = 0; i < J; ++i) {
              std::string jackknifeName = std::string(name) + "_jackknife_partial_" + std::to_string(i);
              std::string jackknifeTitle = std::string(title) + " Jackknife Partial " + std::to_string(i);
              jackknifeHistsVec_.emplace_back(std::make_shared<TH>(jackknifeName.c_str(), jackknifeTitle.c_str(), nbin, xmin, xmax));
            }
          }
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

  // Jacknife Average Hist
  std::shared_ptr<TH> getJackknifeAverage() {
    if (!jackknifeAverage_) {
      throw std::runtime_error("Jackknife average histogram is not available. Please call getCovJackKnife() first.");
    }
    return jackknifeAverage_;
  }
  // Jackknife covariance method
  TH2D* getCovJackKnife(){
    if (!hasCovJackKnife_) {
      hasCovJackKnife_ = true;
      // Create the covariance matrix as a 2D histogram
      jackknifeCov_ = std::make_shared<Result2D_t>(finalHistogram_->GetName(), finalHistogram_->GetTitle(),
                                          finalHistogram_->GetNbinsX(), finalHistogram_->GetXaxis()->GetXmin(),
                                          finalHistogram_->GetXaxis()->GetXmax(),
                                          finalHistogram_->GetNbinsX(), finalHistogram_->GetXaxis()->GetXmin(),
                                          finalHistogram_->GetXaxis()->GetXmax());
      // make all the J jackknife histograms by summing all but one of the jackknifeHistsVec_
      std::vector<std::shared_ptr<TH>> jackknifeHistsVecCopy ;
      jackknifeHistsVecCopy.reserve(jackknifeHistsVec_.size());
      for (int i = 0; i < J; ++i) {
        // Create a temporary histogram to hold the sum of all but one jackknife histogram
        jackknifeHistsVecCopy.push_back(std::make_shared<TH>(*finalHistogram_));
        // reset the last histogram
        auto tempHist = jackknifeHistsVecCopy.back();
        tempHist->Reset();
        for (int j = 0; j < J; ++j) {
          if (j != i) {
            tempHist->Add(jackknifeHistsVec_[j].get());
          }
        }
        // scale the Histogram by the missing data factor J/(J-1)
        tempHist->Scale(static_cast<double>(J) / (J - 1));
      }
      // compute the average of the jackknife histograms and save it in the class member shared ptr
      jackknifeAverage_ = std::make_shared<TH>(*finalHistogram_);
      jackknifeAverage_->Reset();
      for (const auto &hist : jackknifeHistsVecCopy) {
        jackknifeAverage_->Add(hist.get());
      }
      jackknifeAverage_->Scale(1.0 / J);
      auto& averageHist = jackknifeAverage_;
      // TODO: check that averageHist is == to finalHistogram. 

      /*
      std::cout <<"--------DEBUG----------"<<std::endl;
      for(int i=1 ;i<=averageHist->GetNbinsX();++i){
        // print averageHist and finalHistogram
        std::cout << "Average Hist Bin " << i << ": " << averageHist->GetBinContent(i) 
                  << ", Final Hist Bin " << i << ": " << finalHistogram_->GetBinContent(i) << std::endl;
      }
      std::cout <<"------------------"<<std::endl;
      */

      // Fill the covariance matrix with the jackknife resampling histograms as: (J-1) / J * sum (y_j - average)* (y_j-average)^T
      for (const auto &hist : jackknifeHistsVecCopy) {
        // Compute the difference (y_j - average)
        auto diffHist = std::make_shared<TH>(*hist);
        diffHist->Add(averageHist.get(), -1.0);
        
        // Compute the outer product (y_j - average) * (y_j - average)^T
        for (int i = 1; i <= diffHist->GetNbinsX(); ++i) {
          for (int j = 1; j <= diffHist->GetNbinsX(); ++j) {
            double val_i = diffHist->GetBinContent(i);
            double val_j = diffHist->GetBinContent(j);
            double outerProduct = val_i * val_j * (J - 1.0) / J;
            jackknifeCov_->Fill(diffHist->GetBinCenter(i), diffHist->GetBinCenter(j), outerProduct);
          }
        }
      }
    }
    return jackknifeCov_.get();

  }

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
          if (J>1){
            //std::cout<<" FILL FLUSH:" <<eventId << " J: " << J << std::endl; // DEBUG
            hasCovJackKnife_ = false;//here? it should be always reset if I fill something else
            //if all and J>1 throw exception ValueError
            if (all) throw std::runtime_error("OversampledTH: Cannot flush all with jackknife resampling enabled. Please flush individual events.");
            
            // Save histogram for jackknife BEFORE calling fillPartialHists_ (which deletes the histograms)
            auto minJackknifeN = std::min_element(jackknifeN_.begin(), jackknifeN_.end());
            int minJackknifeIndex = std::distance(jackknifeN_.begin(), minJackknifeN);

            //std::cout<<" FILL JACK FLUSH:" <<eventId << " J: " << J << " minJackknifeIndex: " << minJackknifeIndex << std::endl; // DEBUG
            
            // Safety checks and save histogram before it gets deleted
            TH* histToAdd = nullptr;
            if (!it->second.empty()) {
              auto firstSlotHist = it->second.begin()->second;
              if (firstSlotHist && minJackknifeIndex < jackknifeHistsVec_.size() && jackknifeHistsVec_[minJackknifeIndex]) {
                //std::cout << "Adding histogram with " << firstSlotHist->GetEntries() << " entries to jackknife " << minJackknifeIndex << std::endl;
                // Clone the histogram before it gets deleted
                histToAdd = (TH*)firstSlotHist->Clone();
              }
            }
            
            //fillPartialHists_(it->second);
            
            // Add the cloned histogram to jackknife
            if (histToAdd) {
              jackknifeHistsVec_[minJackknifeIndex]->Add(histToAdd);
              jackknifeN_[minJackknifeIndex] += 1;
              delete histToAdd; // Clean up the temporary clone
              //std::cout<<" FILL JACK FLUSH DONE:" <<eventId << " J: " << J << std::endl; // DEBUG
            }
          }
          // Now call fillPartialHists which will delete the originals 
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
