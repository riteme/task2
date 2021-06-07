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

class Index {
public:
    Index();

    auto size() -> int {
        return _n_appended;
    }

    void append(int c);
    void append(const BioSeq &s);
    void build();
    auto rpset(const Token &t) -> std::vector<int>;

    auto next(const Token &t, int c) -> Token;
    auto locate(const BioSeq &s) -> Token;
    auto align(const BioSeq &s) -> Alignment;
    auto fuzzy_locate(const BioSeq &s) -> Location;

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
