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

#include "alias.h"
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
    Analyser(const std::string &sequences_file_path,
             std::uint8_t kmer_len,
             std::uint8_t window_len
    );
    Analyser(const std::string &sequences_file_path,
             std::uint8_t kmer_len,
             std::uint8_t window_len,
             int start,
             int end
    );
    [[nodiscard]] double FindPrecision() const;
    [[nodiscard]] double FindRecall() const;
    SetOverlaps FindTrueRamOverlaps();
    SetOverlaps FindFalsePositive();
    SetOverlaps FindFalseNegative();
    void FindAllTrueOverlaps();
    std::uint32_t num_of_ram_overlaps = 0;
    std::uint32_t num_of_true_ram_overlaps = 0;
    std::uint32_t num_of_true_overlaps = 0;

  private:
    const std::string &path_;
    std::uint8_t kmer_len_;
    std::uint8_t window_len_;
    int start_;
    int end_;
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> targets_;
    SetOverlaps ram_overlaps_;
    SetOverlaps all_true_overlaps_;
    void Initialise();
    bool IsTrueOverlap(std::unique_ptr<biosoup::NucleicAcid> &lhs, std::unique_ptr<biosoup::NucleicAcid> &rhs);
    bool SameContig(std::uint32_t lhs_id, std::uint32_t rhs_id);
  public:
    std::pair<std::uint32_t,
              std::uint32_t> FindRange(std::uint32_t id); // Extract the start and end positions of the simulated read
};

} // namespace ram_analyser

#endif //RAMANALYSER_SRC_ANALYSER_ANALYSER_H_
