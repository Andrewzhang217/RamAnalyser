#include "processor.h"

#include <cstdint>
#include <utility>
#include <iostream>

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace ram_analyser {
    Processor::Processor(std::shared_ptr<thread_pool::ThreadPool> &pool,
                         std::uint8_t kmer_len,
                         std::uint8_t window_len,
                         std::vector<std::unique_ptr<biosoup::NucleicAcid>> &targets) :
            pool_(pool),
            minimizer_engine_(pool, kmer_len, window_len),
            targets_(targets) {
        minimizer_engine_.Minimize(targets_.begin(), targets_.end(), false); // set dafault minhash to false
        minimizer_engine_.Filter(0.001); // set default frequency to 0.001
    }

    std::vector<std::vector<biosoup::Overlap>> Processor::FindAvaOverlaps() {
        return FindOverlaps(0, static_cast<int> (targets_.size()));
    }

    std::vector<std::vector<biosoup::Overlap>> Processor::FindOverlaps(int start, int end) {
        std::uint32_t num_targets = biosoup::NucleicAcid::num_objects;
        biosoup::NucleicAcid::num_objects = 0;
        std::vector<std::vector<biosoup::Overlap>> results;
        for (auto it = start; it < end; ++it) {
            if (it >= num_targets) continue;
            results.emplace_back(minimizer_engine_.Map(targets_[it], true, false, true));
        }
        return results;
    }
} // namespace ram_analyser
