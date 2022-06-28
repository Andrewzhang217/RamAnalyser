#include "analyser.h"
#include "biosoup/progress_bar.hpp"

#include <cstdlib>
#include <utility>

namespace ram_analyser {

Analyser::Analyser(const std::string &sequences_file_path,
                   std::uint8_t kmer_len,
                   std::uint8_t window_len)
        : path(sequences_file_path), kmer_len_(kmer_len), window_len_(window_len) { Initialise(); }
void Analyser::Initialise() {
    std::uint32_t num_threads = 1;
    Input target{path};
    auto targets = target.Sequences();
    auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);
    Processor processor
            {thread_pool, kmer_len_, window_len_, targets};
    overlaps_ = processor.FindOverlaps();
}
Analyser::TrueOverlaps Analyser::FindTrueOverlaps() {
    std::set<std::pair<std::uint32_t, u_int32_t>> set;
    for (const auto &it : overlaps_) {
        for (const auto &overlap : it) {
            if (set.find(std::pair(overlap.lhs_id, overlap.rhs_id)) != set.end()) continue;
            if (set.find(std::pair(overlap.rhs_id, overlap.lhs_id)) != set.end()) continue;
            if (IsTrueOverlap(overlap)) {
                set.emplace(overlap.lhs_id, overlap.rhs_id);
            }
        }
    }
    num_of_true_overlaps = set.size();
    return {std::move(set)};
}
double Analyser::FindPrecision() {
    return static_cast<double> (num_of_true_overlaps) / static_cast<double>(overlaps_.size());
}
double Analyser::FindRecall() {
    // tp / (tp + fn)
    return 0;
}
bool Analyser::IsTrueOverlap(const biosoup::Overlap &overlap) {
    std::uint32_t lhs_id = overlap.lhs_id;
    std::uint32_t rhs_id = overlap.rhs_id;
    // targets do not contain lhs read
    if (targets_.size() <= lhs_id || targets_[lhs_id] == nullptr) {
        return false;
    }
    // targets do not contain rhs read
    if (targets_.size() <= rhs_id || targets_[rhs_id] == nullptr) {
        return false;
    }
    // check whether they really overlap
    std::pair left = FindRange(lhs_id);
    std::pair right = FindRange(rhs_id);
    return (left.second >= right.first && right.second >= left.first)
            || (right.second >= left.first && left.second >= right.first);
}
std::pair<std::uint32_t, std::uint32_t> Analyser::FindRange(uint32_t id) {
    std::string name = targets_[id]->name;
    // todo: extract start and end pos of read from its name string
    return {0, 0};
}
}
