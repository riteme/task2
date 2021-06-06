#include <cstdio>

#include <map>
#include <queue>
#include <tuple>

#include "index.hpp"


namespace {

struct Key {
    int x = 1, y = 0;

    bool operator<(const Key &rhs) const {
        return std::tie(x, y) < std::tie(rhs.x, rhs.y);
    }
};

struct Value {
    int t = 0, l = 0;

    bool operator<(const Value &rhs) const {
        return std::make_tuple(t, -l) < std::make_tuple(rhs.t, -rhs.l);
    }
};

static int current_len;

struct Pair {
    Key key;
    Value value;

    auto estimated() const -> Value {
        // TODO: why (t + 2*r)-heuristic is not effective?

        int remain = current_len - key.y;

        return {
            value.t + 2 * remain,
            value.l
        };

        // return value;
    }

    bool operator<(const Pair &rhs) const {
        return rhs.estimated() < estimated();
        // return rhs.value < value;
    }
};

}

namespace core {

auto Index::align(const BioSeq &s) -> Alignment {
    int n = s.size();
    int max_reach = 0;
    current_len = n;

    size_t max_size = 0;
    std::priority_queue<Pair> q;
    q.push(Pair());

    std::map<Key, Value> f;
    f[Key()] = Value();

    auto probe = [&](const Pair &t) {
        auto it = f.find(t.key);

        if (it == f.end() || t.value < it->second) {
            if (it == f.end())
                f[t.key] = t.value;
            else
                it->second = t.value;

            q.push(t);
        }
    };

    Pair opt;
    do {
        max_size = std::max(max_size, q.size());
        auto u = q.top();
        max_reach = std::max(max_reach, u.value.l + u.key.y);
        if (f.size() % 1000000 <= 2)
            fprintf(stderr, "f.size()=%zu, x=%d, y=%d, t=%d, l=%d, max=%d\n",
                f.size(), u.key.x, u.key.y, u.value.t, u.value.l, max_reach);
        if (f.size() > 20000000)
            break;
        q.pop();

        if (f[u.key] < u.value)
            continue;

        if (u.key.y == n) {
            opt = u;
            break;
        }

        auto [x, y] = u.key;
        auto [t, l] = u.value;

        probe({{x, y + 1}, {t + 11, l}});

        for (int c = 0; c < SIGMA; c++) {
            int z = m[x].transition[c];
            if (!z)
                continue;

            probe({{z, y + 1}, {t + (c == CMAP[s[y + 1]] ? 1 : 6), l + 1}});
            probe({{z, y}, {t + 11, l + 1}});
        }
    } while (!q.empty());

    // 20x + (y + l + x) = 2t
    int loss = (2 * opt.value.t - opt.key.y - opt.value.l) / 21;
    printf("opt: x=%d, y=%d, t=%d, l=%d, loss=%d\n",
        opt.key.x, opt.key.y, opt.value.t, opt.value.l, loss);
    printf("q.max_size=%zu, q.size()=%zu, f.size()=%zu\n",
        max_size, q.size(), f.size());

    return Alignment{
        {opt.key.x, opt.value.l},
        opt.value.t
    };
}

}
