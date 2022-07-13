#include <string>
#include <getopt.h>

#include "analyser.h"

static struct option options[] = {
        {"kmer-length", required_argument, nullptr, 'k'},
        {"window-length", required_argument, nullptr, 'w'},
        {"start", required_argument, nullptr, 's'},
        {"end", required_argument, nullptr, 'e'},
        {nullptr, 0, nullptr, 0}
};

int main(int argc, char **argv) {

    std::uint8_t kmer_len = 15;
    std::uint8_t window_len = 5;
    int start = 0;
    int end = 0;

    const char *optstr = "k:w:s:e";
    char arg;
    while ((arg = static_cast<char>(getopt_long(argc, argv, optstr, options, nullptr))) != -1) {
        switch (arg) {
            case 'k': kmer_len = std::stoi(optarg);
                break;
            case 'w': window_len = std::stoi(optarg);
                break;
            case 's': start = std::stoi(optarg);
                break;
            case 'e': end = std::stoi(optarg);
                break;
            default: return 1;
        }
    }
    std::string sequence_file_path = argv[optind];
    ram_analyser::Analyser analyser{sequence_file_path, kmer_len, window_len, start, end};
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
