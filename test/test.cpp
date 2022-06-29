#include "../src/analyser.h"

#include <algorithm>
#include <iterator>
#include <string_view>

#include "catch2/catch_test_macros.hpp"

namespace ram_analyser_testing {
static const std::string
        &input_file_path = "../test/nanosim_drosophila_ref_iso-profile_100k_extend_aligned_reads.fasta";
static const int k_mer_length = 5;
static const int window_length = 5;
ram_analyser::Analyser analyser{input_file_path, k_mer_length, window_length};
TEST_CASE("Unit test: Input") {
    ram_analyser::Input input{input_file_path};
    auto targets = input.Sequences();
    REQUIRE(!targets.empty());
}
TEST_CASE("Unit test: FindAllTrueOverlaps") {
    analyser.FindAllTrueOverlaps();
}
TEST_CASE("Unit test: FindTrueRamOverlaps") {}
TEST_CASE("Unit test: FindPrecision") {}
TEST_CASE("Unit test: FindRecall") {}
TEST_CASE("Component test") {

}
} // namespace ram_analyser_testing
