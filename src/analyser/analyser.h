#ifndef RAMANALYSER_SRC_ANALYSER_ANALYSER_H_
#define RAMANALYSER_SRC_ANALYSER_ANALYSER_H_

#include <atomic>
#include <memory> // unique_ptr
#include <utility>
#include <vector>
#include <iostream>

#include "align_result.h"
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "../edlib/edlib.h"
#include "../ram/minimizer_engine.h"
#include "thread_pool/thread_pool.hpp"

namespace ram_analyser {
class Analyser {

  public:
    Analyser(std::vector<std::string> sequences_file_paths,
             std::shared_ptr<thread_pool::ThreadPool> &pool,
             std::uint8_t kmer_len,
             std::uint8_t window_len);

    std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(const std::string &path);
    static bool IsSuffix(const std::string &s, const std::string &suff);
    std::vector<biosoup::Overlap> FindOverlaps();

    void find_true_overlaps();

    void find_RAM_overlaps();

    void within_each();

    void true_positives();

    void true_positives_align_part();

    void false_positives();

    void false_positives_align_part();

    void false_negatives();

    void RAM_overlaps_true_reads();

    void RAM_overlaps_simulated_reads();

    void run();

  private:
    std::vector<std::string> paths_;
    ram::MinimizerEngine minimizer_engine_;

};

} // namespace ram_analyser

#endif //RAMANALYSER_SRC_ANALYSER_ANALYSER_H_
