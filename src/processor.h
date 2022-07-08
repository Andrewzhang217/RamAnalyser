#ifndef RAMANALYSER_SRC_ANALYSER_PROCESSOR_H_
#define RAMANALYSER_SRC_ANALYSER_PROCESSOR_H_

#include "alias.h"
#include "ram/minimizer_engine.hpp"

#include <cstdint>

namespace ram_analyser {
class Processor {

  public:
    Processor(std::shared_ptr<thread_pool::ThreadPool> &pool,
              std::uint8_t kmer_len,
              std::uint8_t window_len,
              std::vector<std::unique_ptr<biosoup::NucleicAcid>> &targets);

    SetOverlaps FindAvaOverlaps();
    SetOverlaps FindOverlaps(int start,
                             int end); // find the overlaps among the reads within the range(start, end) and all reads
    static SetOverlaps ConvertRamOverlapsToIds(std::vector<std::vector<biosoup::Overlap>>& overlaps);
  private:
    ram::MinimizerEngine minimizer_engine_;
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> &targets_;
};

}// namespace ram_analyser

#endif //RAMANALYSER_SRC_ANALYSER_PROCESSOR_H_
