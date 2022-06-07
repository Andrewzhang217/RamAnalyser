#include "analyser.h"
#include "biosoup/progress_bar.hpp"

#include <cstdlib>

namespace ram_analyser {
Analyser::Analyser(std::vector<std::string> sequences_file_paths,
                   std::shared_ptr<thread_pool::ThreadPool> &pool,
                   std::uint8_t kmer_len,
                   std::uint8_t window_len,
                   std::uint32_t num_threads,
                   bool minhash,
                   double frequency)
        : paths_(std::move(sequences_file_paths)),
          minimizer_engine_(pool, kmer_len, window_len),
          num_threads_(num_threads),
          minhash_(minhash),
          frequency_(frequency) {
}
void Analyser::find_true_overlaps() {

}
void Analyser::find_RAM_overlaps() {

}
void Analyser::within_each() {

}
void Analyser::true_positives() {

}
void Analyser::true_positives_align_part() {

}
void Analyser::false_positives() {

}
void Analyser::false_positives_align_part() {

}
void Analyser::false_negatives() {

}
void Analyser::RAM_overlaps_true_reads() {

}
void Analyser::RAM_overlaps_simulated_reads() {

}
void Analyser::run() {

}
std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> Analyser::CreateParser(const std::string &path) {

    if (IsSuffix(path, ".fasta") || IsSuffix(path, ".fasta.gz") ||
            IsSuffix(path, ".fna") || IsSuffix(path, ".fna.gz") ||
            IsSuffix(path, ".fa") || IsSuffix(path, ".fa.gz")) {
        try {
            return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path);  // NOLINT
        } catch (const std::invalid_argument &exception) {
            std::cerr << exception.what() << std::endl;
            return nullptr;
        }
    }
    if (IsSuffix(path, ".fastq") || IsSuffix(path, ".fastq.gz") ||
            IsSuffix(path, ".fq") || IsSuffix(path, ".fq.gz")) {
        try {
            return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(path);  // NOLINT
        } catch (const std::invalid_argument &exception) {
            std::cerr << exception.what() << std::endl;
            return nullptr;
        }
    }
    return nullptr;
}
bool Analyser::IsSuffix(const std::string &s, const std::string &suff) {
    return s.size() >= suff.size() && s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
}
std::vector<std::unique_ptr<biosoup::NucleicAcid>> Analyser::FindOverlaps() {

    auto target_parser = CreateParser(paths_[0]);
    bool is_all_versus_all = false;
    std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> sequence_parser;
    if (paths_.size() > 1) {
        sequence_parser = CreateParser(paths_[1]);
        is_all_versus_all = paths_[0] == paths_[1];
    } else {
        sequence_parser = CreateParser(paths_[0]);
        is_all_versus_all = true;
    }

    auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads_);
    biosoup::Timer timer{};
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;

    while (true) {
        timer.Start();

        std::vector<std::unique_ptr<biosoup::NucleicAcid>> targets;
        try {
            targets = target_parser->Parse(1ULL << 32);
        } catch (std::invalid_argument &exception) {
            std::cerr << exception.what() << std::endl;
            return {};
        }

        if (targets.empty()) {
            break;
        }

        std::cerr << "[ram::] parsed " << targets.size() << " targets "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;

        timer.Start();

        minimizer_engine_.Minimize(targets.begin(), targets.end(), minhash_);
        minimizer_engine_.Filter(frequency_);

        std::cerr << "[ram::] minimized targets "
                  << std::fixed << timer.Stop() << "s"
                  << std::endl;

        std::uint64_t num_targets = biosoup::NucleicAcid::num_objects;
        biosoup::NucleicAcid::num_objects = 0;

        while (true) {
            timer.Start();
            try {
                sequences = sequence_parser->Parse(1U << 30);
            } catch (std::invalid_argument &exception) {
                std::cerr << exception.what() << std::endl;
                return {};
            }

            if (sequences.empty()) {
                break;
            }

            std::vector<std::future<std::vector<biosoup::Overlap>>> futures;
            for (const auto &it : sequences) {
                if (is_all_versus_all && it->id >= num_targets) {
                    continue;
                }
                futures.emplace_back(thread_pool->Submit(
                        [&](const std::unique_ptr<biosoup::NucleicAcid> &sequence)
                                -> std::vector<biosoup::Overlap> {
                            return minimizer_engine_.Map(sequence, is_all_versus_all, is_all_versus_all, minhash_);
                        },
                        std::ref(it)));
            }

            biosoup::ProgressBar bar{
                    static_cast<std::uint32_t>(futures.size()), 16};

            std::uint64_t rhs_offset = targets.front()->id;
            std::uint64_t lhs_offset = sequences.front()->id;
            for (auto &it : futures) {
                for (const auto &jt : it.get()) {
                    std::cout << sequences[jt.lhs_id - lhs_offset]->name << "\t"
                              << sequences[jt.lhs_id - lhs_offset]->inflated_len << "\t"
                              << jt.lhs_begin << "\t"
                              << jt.lhs_end << "\t"
                              << (jt.strand ? "+" : "-") << "\t"
                              << targets[jt.rhs_id - rhs_offset]->name << "\t"
                              << targets[jt.rhs_id - rhs_offset]->inflated_len << "\t"
                              << jt.rhs_begin << "\t"
                              << jt.rhs_end << "\t"
                              << jt.score << "\t"
                              << std::max(
                                      jt.lhs_end - jt.lhs_begin,
                                      jt.rhs_end - jt.rhs_begin) << "\t"
                              << 255
                              << std::endl;
                }

                if (++bar) {
                    std::cerr << "[ram::] mapped " << bar.event_counter() << " sequences "
                              << "[" << bar << "] "
                              << std::fixed << timer.Lap() << "s"
                              << "\r";
                }
            }
            std::cerr << std::endl;
            timer.Stop();

            if (is_all_versus_all && biosoup::NucleicAcid::num_objects >= num_targets) {
                break;
            }
        }

        sequence_parser->Reset();
        biosoup::NucleicAcid::num_objects = num_targets;
    }

    std::cerr << "[ram::] " << timer.elapsed_time() << "s" << std::endl;
    return sequences;
}

}
