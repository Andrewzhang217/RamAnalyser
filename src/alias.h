#ifndef RAMANALYSER_ALIAS_H
#define RAMANALYSER_ALIAS_H

#include <unordered_set>

namespace ram_analyser {
struct pred {
    bool operator()(const std::pair<std::uint32_t, std::uint32_t> &pair1,
                    const std::pair<std::uint32_t, std::uint32_t> &pair2) const {
        if (pair1.first == pair2.first) return pair1.second == pair2.second;
        if (pair1.second == pair2.first) return pair1.first == pair2.second;
        return false;
    }
};
struct hash {
    std::size_t operator()(const std::pair<std::uint32_t, std::uint32_t> &pair) const {
        return std::hash<std::uint32_t>{}(pair.first) ^ std::hash<std::uint32_t>{}(pair.second);
    }
};
using SetOverlaps = std::unordered_set<std::pair<std::uint32_t, u_int32_t>, hash, pred>;
} // namespace ram_analyser
#endif //RAMANALYSER_ALIAS_H
