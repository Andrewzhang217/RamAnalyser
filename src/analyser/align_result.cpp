#include "align_result.h"

namespace ram_analyser {

AlignResult::AlignResult(AlignResult &&r) noexcept
        : target_columns(std::move(r.target_columns)), ins_columns(std::move(r.ins_columns)),
          width(r.width) {}
AlignResult &AlignResult::operator=(AlignResult &&r) noexcept {

    if (this != &r) {
        this->target_columns = std::move(r.target_columns);
        this->ins_columns = std::move(r.ins_columns);
        this->width = r.width;
    }
    return *this;

}

}
