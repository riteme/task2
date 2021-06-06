#pragma once

#include <limits>

#include "common.hpp"


namespace core {

constexpr int SIGMA = 5;
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

struct Alignment {
    Token token;
    int loss;
};

class Index {
public:
    Index();

    void append(int c);
    void append(const BioSeq &s);
    void build();
    auto rpset(const Token &t) -> std::vector<int>;

    auto next(const Token &t, int c) -> Token;
    auto locate(const BioSeq &s) -> Token;
    auto align(const BioSeq &s) -> Alignment;

private:
    struct Node {
        int maxlen = 0, fail = 0;
        int index = 0;
        int transition[SIGMA] = {0};

        struct {
            int in = 0, out = 0;
        } dfn;
        std::vector<int> children;
    };

    int _last = 1;
    std::vector<Node> m;
    std::vector<int> sorted;

    auto _allocate(int n) -> int;
    void _copy(int dst, int src);
    auto _append(int x, int c) -> int;
    void _traverse(int x, int &count);
};

}
