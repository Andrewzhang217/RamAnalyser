#ifndef RAMANALYSER_SRC_ANALYSER_PROCESSOR_H_
#define RAMANALYSER_SRC_ANALYSER_PROCESSOR_H_

#include "ram/minimizer_engine.hpp"

#include <cstdint>

namespace ram_analyser {
class Processor {

  public:
    Processor(std::shared_ptr<thread_pool::ThreadPool> &pool,
              std::uint8_t kmer_len,
              std::uint8_t window_len,
              std::vector<std::unique_ptr<biosoup::NucleicAcid>> &targets);

    std::vector<std::vector<biosoup::Overlap>> FindAvaOverlaps();
    std::vector<std::vector<biosoup::Overlap>> FindOverlaps(int start,
                                                            int end); // find the overlaps among the reads within the range(start, end) and all reads

  private:
    std::shared_ptr<thread_pool::ThreadPool> &pool_;
    ram::MinimizerEngine minimizer_engine_;
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &targets_;
};

}// namespace ram_analyser

#endif //RAMANALYSER_SRC_ANALYSER_PROCESSOR_H_
