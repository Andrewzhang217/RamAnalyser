#include "../src/analyser.h"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <string_view>

#include "catch2/catch_test_macros.hpp"

namespace ram_analyser_testing {
static const std::string
        &input_file_path = "../test/nanosim_drosophila_ref_iso-profile_100k_extend_aligned_reads.fasta";
static const int k_mer_length = 5;
static const int window_length = 5;
//ram_analyser::Analyser analyser{input_file_path, k_mer_length, window_length};
TEST_CASE("Unit test: Input") {
    ram_analyser::Input input{input_file_path};
    auto targets = input.Sequences();
    REQUIRE(!targets.empty());
}
TEST_CASE("Unit test: Processor") {
    ram_analyser::Input input{input_file_path};
    auto targets = input.Sequences();
    auto thread_pool = std::make_shared<thread_pool::ThreadPool>(10);
    ram_analyser::Processor processor{thread_pool, 5, 5, targets};
    auto overlaps_first_1k = processor.FindOverlaps(0, 1000);
    std::cout << overlaps_first_1k.size() << std::endl;
    REQUIRE(!overlaps_first_1k.empty());
}
TEST_CASE("Unit test: Analyser") {
    ram_analyser::Analyser analyser{input_file_path, 5, 5, 0, 1000};
    SECTION("Unit test: FindTrueRamOverlaps") {
        auto true_ram_overlaps = analyser.FindTrueRamOverlaps();
        std::cout << analyser.num_of_true_ram_overlaps << std::endl;
        if (analyser.num_of_true_ram_overlaps > 0) {
            REQUIRE(true_ram_overlaps.size() == analyser.num_of_true_ram_overlaps);
        }
    }SECTION("Unit test: FindPrecision") {

    }SECTION("Unit test: FindRecall") {

    }
}
TEST_CASE("Component test") {

}
} // namespace ram_analyser_testing
