#ifndef RAMANALYSER_SRC_ANALYSER_INPUT_H_
#define RAMANALYSER_SRC_ANALYSER_INPUT_H_

#include <vector>
#include <memory>
#include "biosoup/nucleic_acid.hpp"

namespace ram_analyser {
class Input {

  public:
    explicit Input(const std::string &path);
    static bool IsSuffix(const std::string &s, const std::string &suff);
    std::vector<std::unique_ptr<biosoup::NucleicAcid>> ParseFile();

  private:
    const std::string &path_;
};
}

#endif //RAMANALYSER_SRC_ANALYSER_INPUT_H_
