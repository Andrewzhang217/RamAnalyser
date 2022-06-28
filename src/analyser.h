#ifndef RAMANALYSER_SRC_ANALYSER_ANALYSER_H_
#define RAMANALYSER_SRC_ANALYSER_ANALYSER_H_

#include <atomic>
#include <cstdint>
#include <memory> // unique_ptr
#include <map>
#include <utility>
#include <set>
#include <vector>
#include <iostream>

#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "input.h"
#include "processor.h"
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

namespace ram_analyser {
class Analyser {
  public:
    using TrueOverlaps = std::set<std::pair<std::uint32_t, u_int32_t>>;
    Analyser(const std::string &sequences_file_path,
             std::uint8_t kmer_len,
             std::uint8_t window_len
    );

    double FindPrecision();
    double FindRecall();
    TrueOverlaps FindTrueOverlaps();
    static std::uint32_t num_of_true_overlaps;

  private:
    const std::string &path;
    std::uint8_t kmer_len_;
    std::uint8_t window_len_;
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> targets_;
    std::vector<std::vector<biosoup::Overlap>> overlaps_;
    void Initialise();
    bool IsTrueOverlap(const biosoup::Overlap &overlap);
    std::pair<std::uint32_t,
              std::uint32_t> FindRange(uint32_t id); // Extract the start and end positions of the simulated read
};

} // namespace ram_analyser

#endif //RAMANALYSER_SRC_ANALYSER_ANALYSER_H_