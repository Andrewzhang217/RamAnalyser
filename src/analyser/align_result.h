#ifndef RAMANALYSER_SRC_ANALYSER_ALIGN_RESULT_H_
#define RAMANALYSER_SRC_ANALYSER_ALIGN_RESULT_H_

#include <cstdint>
#include <vector>

namespace ram_analyser {
class AlignResult {

  public :
    AlignResult() = default;
    AlignResult(AlignResult &&r) noexcept;
    AlignResult &operator=(AlignResult &&r) noexcept;

  private:
    std::vector<std::vector<char>> target_columns;
    // which target position? -> which ins at that position? -> which row?
    std::vector<std::vector<std::vector<char>>>
            ins_columns;
    std::uint32_t width{};

};
}

#endif //RAMANALYSER_SRC_ANALYSER_ALIGN_RESULT_H_
