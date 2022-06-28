#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "bioparser/parser.hpp"
#include "input.h"

#include <iostream>
#include <memory>

namespace ram_analyser {
Input::Input(const std::string &path) :
        path_(path) {
    ParseFile();
}

void Input::ParseFile() {
    if (IsSuffix(path_, ".fasta") || IsSuffix(path_, ".fasta.gz") ||
            IsSuffix(path_, ".fna") || IsSuffix(path_, ".fna.gz") ||
            IsSuffix(path_, ".fa") || IsSuffix(path_, ".fa.gz")) {
        try {
            auto parser = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path_);
            sequences_ = parser->Parse(-1);
        } catch (const std::invalid_argument &exception) {
            std::cerr << exception.what() << std::endl;
            return;
        }
    }
    if (IsSuffix(path_, ".fastq") || IsSuffix(path_, ".fastq.gz") ||
            IsSuffix(path_, ".fq") || IsSuffix(path_, ".fq.gz")) {
        try {
            auto parser = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(path_);
            sequences_ = parser->Parse(-1);
        } catch (const std::invalid_argument &exception) {
            std::cerr << exception.what() << std::endl;
            return;
        }
    }
    throw std::invalid_argument("Error: Invalid sequences file format.");
}

bool Input::IsSuffix(const std::string &s, const std::string &suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}
std::vector<std::unique_ptr<biosoup::NucleicAcid>> Input::Sequences() {
    return std::move(sequences_);
}
}
