#include <string>
#include <getopt.h>

#include "analyser.h"

void Help() {
    std::cout << "usage: ./RamAnalyser [options ...] <sequences>\n"
                 "\n"
                 "  <target>/<sequences>\n"
                 "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
                 "\n"
                 "  options:\n"
                 "    -k <std::uint8_t>\n"
                 "      default: 15\n"
                 "      length of minimizers\n"
                 "    -w <std::uint8_t>\n"
                 "      default: 5\n"
                 "      length of sliding window from which minimizers are sampled\n"
                 "    -f <double>\n"
                 "      default: 0.001\n"
                 "      threshold for ignoring most frequent minimizers\n"
                 "    -M <bool>\n"
                 "      default: false\n"
                 "      set minhash to true, use only a portion of all minimizers\n"
                 "    -t <std::uint8_t>\n"
                 "      default: 10\n"
                 "      number of threads\n"
                 "    -h\n"
                 "      prints the usage"
                 << std::endl;
}
int main(int argc, char **argv) {

    std::uint8_t kmer_len = 15;
    std::uint8_t window_len = 5;
    int size = 1000;
    bool minhash = false;
    std::uint8_t num_threads = 10;
    double frequency = 0.001;

    const char *optstr = "k:w:s:t:f:Mh";
    char arg;
    while ((arg = static_cast<char>(getopt(argc, argv, optstr))) != -1) {
        switch (arg) {
            case 'k': kmer_len = std::stoi(optarg);
                break;
            case 'w': window_len = std::stoi(optarg);
                break;
            case 's': size = std::stoi(optarg);
                break;
            case 'M': minhash = true;
                break;
            case 't': num_threads = std::stoi(optarg);
                break;
            case 'f': frequency = std::stod(optarg);
                break;
            case 'h': Help(); return 0;
            default: return 1;
        }
    }
    std::string sequence_file_path = argv[optind];
    ram_analyser::Analyser analyser{sequence_file_path, kmer_len, window_len, size, minhash, num_threads, frequency};
    analyser.FindTrueRamOverlaps();
    std::cout << "Number of true positives: " << analyser.num_of_true_ram_overlaps << "\n";
    analyser.FindFalsePositive();
    std::cout << "Number of false positives: " << analyser.num_of_false_positives << "\n";
    analyser.FindFalseNegative();
    std::cout << "Number of false negatives: " << analyser.num_of_false_negatives << "\n";
    std::cout << "Precision: " << analyser.FindPrecision() << std::endl;
    std::cout << "Recall: " << analyser.FindRecall() << std::endl;
    return 0;
}
