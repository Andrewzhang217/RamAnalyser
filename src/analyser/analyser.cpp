#include "analyser.h"

namespace ram_analyser {
Analyser::Analyser(const char **sequences_file_paths,
                   std::shared_ptr<thread_pool::ThreadPool> &pool,
                   std::uint8_t kmer_len,
                   std::uint8_t window_len,
                   double freq) {

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

}
