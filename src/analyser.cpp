#include "analyser.h"
#include "biosoup/progress_bar.hpp"

#include <cstdlib>
#include <utility>

std::uint32_t ram_analyser::Analyser::num_of_true_ram_overlaps{0};
std::uint32_t ram_analyser::Analyser::num_of_true_overlaps{0};

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
Analyser::SetOverlaps Analyser::FindTrueRamOverlaps() {
    FindAllTrueOverlaps();
    std::set<std::pair<std::uint32_t, u_int32_t>> set;
    for (const auto &it : overlaps_) {
        for (const auto &overlap : it) {
            if (all_true_overlaps_.find(std::pair(overlap.lhs_id, overlap.rhs_id)) != all_true_overlaps_.end()) {
                set.emplace(overlap.lhs_id, overlap.rhs_id);
            }
        }
    }
    num_of_true_ram_overlaps = set.size();
    return {std::move(set)};
}
Analyser::SetOverlaps Analyser::FindFalsePositive() {
    return {};
}
void Analyser::FindAllTrueOverlaps() {
    std::set<std::pair<std::uint32_t, u_int32_t>> set;
    for (auto lhs = targets_.begin(); lhs != targets_.end(); ++lhs) {
        for (auto rhs = lhs + 1; rhs != targets_.end(); ++rhs) {
            if (IsTrueOverlap(*lhs, *rhs)) set.emplace(lhs->get()->id, rhs->get()->id);
        }
    }
    num_of_true_overlaps = set.size();
    all_true_overlaps_ = std::move(set);
}
double Analyser::FindPrecision() {
    return static_cast<double> (num_of_true_ram_overlaps) / static_cast<double>(overlaps_.size());
}
double Analyser::FindRecall() {
    // tp / (tp + fn)
    // fn is the number of overlaps in true but not found by RAM
    if (num_of_true_overlaps == 0) return 0;
    return static_cast<double> (num_of_true_ram_overlaps) / static_cast<double> (num_of_true_overlaps);
}
bool Analyser::IsTrueOverlap(std::unique_ptr<biosoup::NucleicAcid> &lhs, std::unique_ptr<biosoup::NucleicAcid> &rhs) {
    std::pair left = FindRange(lhs->id);
    std::pair right = FindRange(rhs->id);
    return (left.second >= right.first && right.second >= left.first)
            || (right.second >= left.first && left.second >= right.first);
}
std::pair<std::uint32_t, std::uint32_t> Analyser::FindRange(uint32_t id) {
    std::string name = targets_[id]->name;
    std::uint32_t start = 0;
    std::uint32_t end = 0;
    auto curr = std::find_if(name.begin(), name.end(), [](auto c) { return c == '_'; });
    for (auto it = ++curr; !std::isdigit(*it); ++it) {
        start *= 10;
        start += *it;
    }
    int underscore_count = 0;
    for (; underscore_count < 5; ++curr) {
        if (*curr == '_') underscore_count++;
    }
    for (; !std::isdigit(*curr); ++curr) {
        end *= 10;
        end += *curr;
    }
    return {start, end};
}
} // namespace ram_analyser
