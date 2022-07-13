#include "analyser.h"
#include "biosoup/progress_bar.hpp"

#include <cassert>
#include <utility>

namespace ram_analyser {
Analyser::Analyser(const std::string &sequences_file_path, std::uint8_t kmer_len, std::uint8_t window_len, int size,
                   bool minhash, std::uint32_t num_threads, double frequency) :
        path_(sequences_file_path), kmer_len_(kmer_len), window_len_(window_len), size_(size){
    Initialise(minhash, num_threads, size, frequency);
}
void Analyser::Initialise(bool minhash, std::uint32_t num_threads, int size, double frequency) {
    Input target{path_};
    targets_ = std::move(target.Sequences());
    auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);
    Processor processor
            {thread_pool, kmer_len_, window_len_, targets_, minhash, frequency};
    ram_overlaps_ = processor.FindOverlaps(size);
    num_of_ram_overlaps = ram_overlaps_.size();
    FindAllTrueOverlaps();
}
SetOverlaps Analyser::FindTrueRamOverlaps() {
    assert(!ram_overlaps_.empty());
    if (all_true_overlaps_.empty()) return {};
    SetOverlaps set;

    for (const auto &overlap : ram_overlaps_) {
        if (all_true_overlaps_.find(std::pair(overlap.first, overlap.second)) != all_true_overlaps_.end()) {
            set.emplace(overlap.first, overlap.second);
        }
    }
    num_of_true_ram_overlaps = set.size();
    return set;
}
SetOverlaps Analyser::FindFalsePositive() {
    assert(!ram_overlaps_.empty());
    if (all_true_overlaps_.empty()) return ram_overlaps_;
    auto false_positive = ram_overlaps_;

    for (auto &overlap : ram_overlaps_) {
        if (all_true_overlaps_.find(overlap) != all_true_overlaps_.end()) {
            false_positive.erase(overlap);
        }
    }
    num_of_false_positives = false_positive.size();
    return false_positive;
}
SetOverlaps Analyser::FindFalseNegative() {
    if (all_true_overlaps_.empty()) return {};
    auto false_negative = all_true_overlaps_;
    auto false_positive = FindTrueRamOverlaps();

    for (auto &overlap : false_positive) {
        if (false_negative.find(overlap) != false_negative.end()) {
            false_negative.erase(overlap);
        }
    }
    num_of_false_negatives = false_negative.size();
    return false_negative;
}
void Analyser::FindAllTrueOverlaps() {
    SetOverlaps set;
    auto end_pos = size_ == 0 ? targets_.end() : targets_.begin() + size_;
    for (auto lhs = targets_.begin(); lhs != end_pos; ++lhs) {
        for (auto rhs = lhs + 1; rhs != targets_.end(); ++rhs) {
            if (IsTrueOverlap(*lhs, *rhs)) set.emplace(lhs->get()->id, rhs->get()->id);
        }
    }
    num_of_true_overlaps = set.size();
    all_true_overlaps_ = std::move(set);
}
double Analyser::FindPrecision() const {
    assert(!ram_overlaps_.empty());
    return static_cast<double> (num_of_true_ram_overlaps) / static_cast<double>(num_of_ram_overlaps);
}
double Analyser::FindRecall() const {
    // tp / (tp + fn)
    // fn is the number of overlaps in true but not found by RAM
    if (num_of_true_overlaps == 0) return 0;
    return static_cast<double> (num_of_true_ram_overlaps) / static_cast<double> (num_of_true_overlaps);
}
bool Analyser::IsTrueOverlap(std::unique_ptr<biosoup::NucleicAcid> &lhs, std::unique_ptr<biosoup::NucleicAcid> &rhs) {
    if (!SameContig(lhs->id, rhs->id)) return false;
    std::pair left = FindRange(lhs->id);
    std::pair right = FindRange(rhs->id);
    return (left.second >= right.first && right.second >= left.first);
}
std::pair<std::uint32_t, std::uint32_t> Analyser::FindRange(std::uint32_t id) {
    std::string name = targets_[id]->name;
    std::uint32_t start = 0;
    std::uint32_t len = 0;
    auto curr = std::find_if(name.begin(), name.end(), [](auto c) { return c == '_'; });
    curr++;
    for (; std::isdigit(*curr); ++curr) {
        start *= 10;
        start += *curr - '0';
    }
    int underscore_count = 1;
    for (; underscore_count < 5; ++curr) {
        if (*curr == '_') underscore_count++;
    }
    for (; std::isdigit(*curr); ++curr) {
        len *= 10;
        len += *curr - '0';
    }
    return {start, start + len};
}

bool Analyser::SameContig(std::uint32_t lhs_id, std::uint32_t rhs_id) {
    std::string name_lhs = targets_[lhs_id]->name;
    std::string name_rhs = targets_[rhs_id]->name;
    auto curr_lhs = std::find_if(name_lhs.begin(), name_lhs.end(), [](auto c) { return c == '_'; });
    auto curr_rhs = std::find_if(name_rhs.begin(), name_rhs.end(), [](auto c) { return c == '_'; });
    return name_lhs.substr(0, curr_lhs - name_lhs.begin()) == name_rhs.substr(0, curr_rhs - name_rhs.begin());
}
} // namespace ram_analyser
