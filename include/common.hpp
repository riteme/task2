#pragma once

#include <cstdint>

#include <vector>
#include <string>


namespace core {

using i8 = std::int8_t;
using u8 = std::uint8_t;
using i64 = std::int64_t;
using u64 = std::uint64_t;

// Sequence is 1-indexed.
template <typename T, typename TContainer>
class Sequence {
public:
    TContainer &internal;

    Sequence(TContainer &_internal)
        : internal(_internal) {}

    auto begin() {
        return internal.begin();
    }

    auto begin() const {
        return internal.begin();
    }

    auto end() {
        return internal.begin();
    }

    auto end() const {
        return internal.end();
    }

    auto size() const -> int {
        return internal.size();
    }

    auto &operator[](int i) {
        return internal[i - 1];
    }

    auto &operator[](int i) const {
        return internal[i - 1];
    }
};

using BioSeq = Sequence<char, std::basic_string<char>>;

}
