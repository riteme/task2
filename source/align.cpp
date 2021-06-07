#include <cstdio>

#include <queue>
#include <tuple>

#include "tsl/robin_map.h"

#include "index.hpp"


namespace {

struct Key {
    int x = 1, y = 0;

    bool operator==(const Key &rhs) const {
        return x == rhs.x && y == rhs.y;
    }

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
            value.t + 3 * remain,
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

namespace std {

template <>
struct hash<Key> {
    auto operator()(const Key &z) const -> size_t {
        return (size_t(z.x) * 19260817) ^ (size_t(z.y) * 0x19260817);
    }
};

};

namespace core {

auto Index::align(const BioSeq &s) -> Alignment {
    int n = s.size();
    current_len = n;

    size_t max_queue_size = 0;
    std::priority_queue<Pair> q;
    q.push(Pair());

    // std::unordered_map<Key, Value> f;
    tsl::robin_map<Key, Value> f;
    f[Key()] = Value();

    auto probe = [&](const Pair &t) {
        auto it = f.find(t.key);

        if (it == f.end() || t.value < it->second) {
            if (it == f.end())
                f[t.key] = t.value;
            else {
                // it->second = t.value;
                it.value() = t.value;
            }

            q.push(t);
        }
    };

    Pair opt;
    do {
        max_queue_size = std::max(max_queue_size, q.size());
        auto u = q.top();
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

            probe({{z, y + 1}, {t + (c == CMAP[s[y + 1]] ? 1 : 11), l + 1}});
            probe({{z, y}, {t + 11, l + 1}});
        }
    } while (!q.empty());

    AlignmentDebugInfo debug;
    debug.n_state_visited = f.size();
    debug.max_queue_size = max_queue_size;

    // 20x + (y + l + x) = 2t
    debug.pure_loss = (2 * opt.value.t - opt.key.y - opt.value.l) / 21;

    return Alignment{
        {opt.key.x, opt.value.l},
        opt.value.t, debug
    };
}

}
