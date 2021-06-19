#pragma once

#include <limits>

#include "common.hpp"


namespace core {

constexpr int ALPHABET_SIZE = 5;
constexpr int CMAP[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

struct Alignment {
    Range range1, range2;
    int loss;
    bool mark;

    auto length1() const -> int {
        return range1.length();
    }

    auto length2() const -> int {
        return range2.length();
    }

    auto rate1() const -> double {
        return double(loss) / length1();
    }

    auto rate2() const -> double {
        return double(loss) / length2();
    }
};

auto full_align(const BioSeq &s1, const BioSeq &s2) -> int;
auto local_align(const BioSeq &s1, const BioSeq &s2) -> Alignment;
auto concat_align(const BioSeq &s1, const BioSeq &s2) -> Alignment;

auto sublocal_span(const BioSeq &s1, const BioSeq &s2) -> Alignment;
auto prefix_span(const BioSeq &s1, const BioSeq &s2) -> Alignment;
auto suffix_span(const BioSeq &s1, const BioSeq &s2) -> Alignment;

class Index {
public:
    struct Token {
        int id = 1;
        int len = 0;
    };

    struct Location {
        bool reversed;
        int left;
        int right;
    };

    struct AlignmentDebugInfo {
        int n_state_visited;
        int max_queue_size;
    };

    struct Alignment {
        Token token;
        int loss;
        AlignmentDebugInfo debug;
    };

    Index();

    auto size() const -> int {
        return _n_appended;
    }

    void append(int c);
    void append(const BioSeq &s);
    void build();

    auto rpset(const Token &t) const -> std::vector<int>;

    auto next(const Token &t, int c) const -> Token;
    auto locate(const BioSeq &s) const -> Token;
    auto align(const BioSeq &s) const -> Alignment;
    auto fuzzy_locate(const BioSeq &s) const -> Location;

private:
    struct Node {
        int maxlen = 0, fail = 0;
        int index = 0;
        int transition[ALPHABET_SIZE] = {0};

        struct {
            int in = 0, out = 0;
        } dfn;
        std::vector<int> children;
    };

    int _last = 1;
    int _n_appended = 0;
    std::vector<Node> m;
    std::vector<int> _sorted;

    auto _allocate(int n) -> int;
    void _copy(int dst, int src);
    auto _append(int x, int c) -> int;
    void _traverse(int x, int &count);
};

}
