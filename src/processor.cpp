#include "processor.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace ram_analyser {
Processor::Processor(std::shared_ptr<thread_pool::ThreadPool> &pool, std::uint8_t kmer_len, std::uint8_t window_len,
                     std::vector<std::unique_ptr<biosoup::NucleicAcid>> &targets, bool minhash, double frequency) :
        minimizer_engine_(pool, kmer_len, window_len),
        targets_(targets),
        minhash_(minhash){
    minimizer_engine_.Minimize(targets_.begin(), targets_.end(), minhash);
    minimizer_engine_.Filter(frequency);
}
SetOverlaps Processor::FindOverlaps(int size) {
    std::uint32_t num_targets = biosoup::NucleicAcid::num_objects;
    biosoup::NucleicAcid::num_objects = 0;
    std::vector<std::vector<biosoup::Overlap>> results;
    for (auto it = 0; it < size; ++it) {
        if (it >= num_targets) continue;
        results.emplace_back(minimizer_engine_.Map(targets_[it], true, false, minhash_));
    }
    return ConvertRamOverlapsToIds(results);
}
SetOverlaps Processor::ConvertRamOverlapsToIds(std::vector<std::vector<biosoup::Overlap>> &overlaps) {
    SetOverlaps set;
    for (const auto &it : overlaps) {
        for (const auto &overlap : it) {
            set.emplace(overlap.lhs_id, overlap.rhs_id);
        }
    }
    return set;
}
} // namespace ram_analyser
