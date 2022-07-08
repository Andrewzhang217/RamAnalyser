#include "../src/analyser.h"

#include <algorithm>
#include <iterator>
#include <iostream>

#include "catch2/catch_test_macros.hpp"

namespace ram_analyser_testing {
static const std::string
        &input_file_path = "../../nanosim_drosophila_ref_iso-profile_100k_extend_aligned_reads.fasta";
static const int k_mer_length = 15;
static const int window_length = 5;
TEST_CASE("Unit test: Input") {
    ram_analyser::Input input{input_file_path};
    auto targets = input.Sequences();
    REQUIRE(!targets.empty());
}

TEST_CASE("Unit test: Minimizer_engine indexing") {
    ram_analyser::Input input{input_file_path};
    auto targets = input.Sequences();
    std::cout << targets.size() << std::endl;
    auto thread_pool = std::make_shared<thread_pool::ThreadPool>(10);
    ram::MinimizerEngine minimizerEngine{thread_pool};
    minimizerEngine.Minimize(targets.begin(), targets.end());
}
TEST_CASE("Unit test: Processor") {
    ram_analyser::Input input{input_file_path};
    auto targets = input.Sequences();
    auto thread_pool = std::make_shared<thread_pool::ThreadPool>(20);
    ram_analyser::Processor processor{thread_pool, 15, 5, targets};
    auto overlaps_first_1k = processor.FindOverlaps(0, 1000);
    std::cout << "num of overlaps found by ram = " << overlaps_first_1k.size() << std::endl;
    REQUIRE(!overlaps_first_1k.empty());
}
TEST_CASE("Unit test: Analyser") {
    ram_analyser::Analyser analyser1{input_file_path, k_mer_length, window_length, 0, 100};
//    SECTION("Unit test: FindAllTrueOverlaps") {
    std::cout << "num of all true overlaps: " << analyser1.num_of_true_overlaps << std::endl;
    //}
    //SECTION("Unit test: FindTrueRamOverlaps") {
    auto true_ram_overlaps = analyser1.FindTrueRamOverlaps();
    std::cout << "num_of_true_ram_overlaps: " << analyser1.num_of_true_ram_overlaps << std::endl; // 3849
    if (analyser1.num_of_true_ram_overlaps > 0) {
        REQUIRE(true_ram_overlaps.size() == analyser1.num_of_true_ram_overlaps);
    }
    auto false_positive = analyser1.FindFalsePositive().size();
    std::cout << "num_of_false_positive_overlaps: " << false_positive << std::endl;
    // }
    //SECTION("Unit test: FindPrecision") {
    std::cout << "num_of_true_ram_overlaps: " << analyser1.num_of_true_ram_overlaps << "\n";
    std::cout << "num_of_ram_overlaps: " << analyser1.num_of_ram_overlaps << "\n";
    auto precision = analyser1.FindPrecision();
    std::cout << "precision: " << precision << std::endl;
    //}
    //SECTION("Unit test: FindRecall") {
    std::cout << "num_of_true_ram_overlaps: " << analyser1.num_of_true_ram_overlaps << "\n";
    std::cout << "num_of_true_overlaps: " << analyser1.num_of_true_overlaps << "\n";
    auto recall = analyser1.FindRecall();
    std::cout << "recall: " << recall << std::endl;
    //}
}
} // namespace ram_analyser_testing
