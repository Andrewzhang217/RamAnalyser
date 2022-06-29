#include "processor.h"

#include <cstdint>
#include <utility>

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace ram_analyser {
Processor::Processor(std::shared_ptr<thread_pool::ThreadPool> &pool,
                     std::uint8_t kmer_len,
                     std::uint8_t window_len,
                     std::vector<std::unique_ptr<biosoup::NucleicAcid>> &targets) :
        pool_(pool),
        minimizer_engine_(pool, kmer_len, window_len),
        targets_(targets) {}
std::vector<std::vector<biosoup::Overlap>> Processor::FindOverlaps() {
    minimizer_engine_.Minimize(targets_.begin(), targets_.end(), false); // set dafault minhash to false
    minimizer_engine_.Filter(0.001); // set default frequency to 0.001
    std::uint32_t num_targets = biosoup::NucleicAcid::num_objects.load();
    biosoup::NucleicAcid::num_objects = 0;
    std::vector<std::vector<biosoup::Overlap>> results;
    while (true) {
        std::vector<std::future<std::vector<biosoup::Overlap>>> futures;
        for (const auto &it : targets_) {
            if (it->id >= num_targets) {
                continue;
            }
            futures.emplace_back(pool_->Submit(
                    [&](const std::unique_ptr<biosoup::NucleicAcid> &sequence)
                            -> std::vector<biosoup::Overlap> {
                        return minimizer_engine_.Map(sequence, true, false, false); // is_ava = true
                    },
                    std::ref(it)));
        }
        if (biosoup::NucleicAcid::num_objects >= num_targets) {
            for (auto &it : futures) {
                results.emplace_back(it.get());
            }
            break;
        }
    }
    return results;
}

} // namespace ram_analyser
