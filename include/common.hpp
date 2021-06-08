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

    Slice(TContainer &_internal)
        : internal(_internal), m_begin(_internal.begin()), m_end(_internal.end()) {}
    Slice(TContainer &_internal, int begin, int end)
        : internal(_internal),
          m_begin(_internal.begin() + (begin - 1)),
          m_end(_internal.begin() + (end - 1)) {}

    auto begin() {
        return m_begin;
    }

    auto begin() const {
        return m_begin;
    }

    auto end() {
        return m_end;
    }

    auto end() const {
        return m_end;
    }

    auto size() const -> int {
        return m_end - m_begin;
    }

    auto &operator[](int i) {
        return *(m_begin + (i - 1));
    }

    auto &operator[](int i) const {
        return *(m_begin + (i - 1));
    }

    auto take(int begin, int end) const -> Slice {
        return {internal, m_begin + (begin - 1), m_begin + (end - 1)};
    }

private:
    using Iterator = decltype(internal.begin());
    Iterator m_begin, m_end;

    Slice(TContainer &_internal, Iterator _begin, Iterator _end)
        : internal(_internal), m_begin(_begin), m_end(_end) {}
};

using BioSeq = Slice<char, std::basic_string<char>>;

}
