#include <algorithm>

#include "index.hpp"


namespace core {

Index::Index() {
    _allocate(2);  // 0 & 1
    m[0].maxlen = -1;
    std::fill(
        std::begin(m[0].transition),
        std::end(m[0].transition),
        1
    );
}

auto Index::_allocate(int n) -> int {
    m.resize(m.size() + n);
    return m.size() - 1;
}

void Index::_copy(int dst, int src) {
    m[dst].fail = m[src].fail;
    m[dst].index = m[src].index;
    std::copy(
        std::begin(m[src].transition),
        std::end(m[src].transition),
        std::begin(m[dst].transition)
    );
}

auto Index::_append(int x, int c) -> int {
    _n_appended++;

    int y = _allocate(1);
    m[y].index = m[y].maxlen = m[x].maxlen + 1;

    while (!m[x].transition[c]) {
        m[x].transition[c] = y;
        x = m[x].fail;
    }

    int p = m[x].transition[c];
    if (m[x].maxlen + 1 != m[p].maxlen) {
        int q = _allocate(1);
        _copy(q, p);
        m[q].maxlen = m[x].maxlen + 1;
        m[p].fail = m[y].fail = q;

        while (m[x].transition[c] == p) {
            m[x].transition[c] = q;
            x = m[x].fail;
        }
    } else
        m[y].fail = p;

    return y;
}

void Index::append(int c) {
    _last = _append(_last, c);
}

void Index::append(const BioSeq &s) {
    for (auto c : s) {
        append(CMAP[static_cast<u8>(c)]);
    }
}

void Index::build() {
    for (int i = 1; i < m.size(); i++) {
        m[i].children.clear();
    }

    for (int i = 2; i < m.size(); i++) {
        int p = m[i].fail;
        m[p].children.push_back(i);
    }

    _sorted.clear();
    _sorted.resize(m.size());
    int count = 0;
    _traverse(1, count);
}

void Index::_traverse(int x, int &count) {
    count++;
    m[x].dfn.in = count;
    _sorted[count] = m[x].index;

    for (int v : m[x].children) {
        _traverse(v, count);
    }

    m[x].dfn.out = count;
}

auto Index::rpset(const Token &t) -> std::vector<int> {
    int l = m[t.id].dfn.in;
    int r = m[t.id].dfn.out + 1;
    std::vector<int> set(_sorted.begin() + l, _sorted.begin() + r);

    std::sort(set.begin(), set.end());
    set.erase(unique(set.begin(), set.end()), set.end());

    return set;
}

auto Index::next(const Token &t, int c) -> Token {
    auto [x, l] = t;
    while (!m[x].transition[c]) {
        x = m[x].fail;
        l = m[x].maxlen;
    }

    return {m[x].transition[c], l + 1};
}

auto Index::locate(const BioSeq &s) -> Token {
    Token t;
    for (int i = 1; i <= s.size(); i++) {
        t = next(t, CMAP[s[i]]);
    }

    return t;
}

}
