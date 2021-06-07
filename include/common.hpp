#pragma once

#include <cstdint>

#include <vector>
#include <string>


namespace core {

using i8 = std::int8_t;
using u8 = std::uint8_t;
using i64 = std::int64_t;
using u64 = std::uint64_t;

bool startswith(const std::string &target, const std::string &pattern);
auto watson_crick_complement(const std::string &s) -> std::string;

// Sequence is 1-indexed.
template <typename T, typename TContainer>
class Slice {
public:
    TContainer &internal;
    int start, length;

    Slice(TContainer &_internal)
        : internal(_internal), start(0), length(_internal.size()) {}
    Slice(TContainer &_internal, int _start, int _length)
        : internal(_internal), start(_start), length(_length) {}

    auto begin() {
        return internal.begin() + start;
    }

    auto begin() const {
        return internal.begin() + start;
    }

    auto end() {
        return internal.begin() + start + length;
    }

    auto end() const {
        return internal.begin() + start + length;
    }

    auto size() const -> int {
        return length;
    }

    auto &operator[](int i) {
        return internal[start + i - 1];
    }

    auto &operator[](int i) const {
        return internal[start + i - 1];
    }

    auto take(int _begin, int _end) const -> Slice {
        return Slice(internal, start + _begin - 1, _end - _begin);
    }
};

using BioSeq = Slice<char, std::basic_string<char>>;

}
