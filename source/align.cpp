#include <cstdio>

#include <queue>
#include <tuple>

#include "tsl/robin_map.h"

#include "index.hpp"


namespace {

constexpr int MISS_COST = 10;
constexpr int CHAR_COST = 1;
constexpr int FULL_COST = MISS_COST + CHAR_COST;
constexpr int H_VALUE = 5;

struct Key {
    int x = 1, y = 0;

    bool operator==(const Key &rhs) const {
        return x == rhs.x && y == rhs.y;
    }

    bool operator<(const Key &rhs) const {
        return std::tie(x, y) < std::tie(rhs.x, rhs.y);
    }
};

struct State {
    Key key;
    int t = 0, l = 0;
};

struct Heuristic {
    int n;

    auto estimate(const State &s) -> int {
        return s.t + H_VALUE * (n - s.key.y);
    }

    bool operator()(const State &lhs, const State &rhs) {
        return estimate(lhs) > estimate(rhs);
    }
};

}

namespace std {

template <>
struct hash<Key> {
    auto operator()(const Key &z) const -> size_t {
        return ((z.x << 15) ^ z.y) * 0x19260817;
    }
};

};

namespace core {

auto local_align(const BioSeq &s1, const BioSeq &s2) -> Alignment {
    constexpr int INF = 0x3f3f3f3f;
    struct Value {
        int t, l;

        static auto zero() -> Value {
            return {0, 0};
        }

        static auto max() -> Value {
            return {INF, INF};
        }

        bool operator<(const Value &rhs) const {
            return t == rhs.t ? l < rhs.l : t < rhs.t;
        }

        auto operator+(const Value &rhs) const -> Value {
            return {t + rhs.t, l + rhs.l};
        }
    };

    auto update = [](Value &dest, const Value &val) {
        if (val < dest)
            dest = val;
    };

    int n = s1.size(), m = s2.size();
    std::vector<Value> f;

    f.resize(m + 1);
    for (int i = 0; i <= m; i++) {
        f[i] = {i, 0};
    }

    auto opt = Value::max();
    int opt_i = 0;
    for (int i = 1; i <= n; i++) {
        for (int j = m; j > 0; j--) {
            f[j] = f[j] + Value{1, 1};
            if (s1[i] == s2[j])
                update(f[j], f[j - 1] + Value{0, 1});
        }

        update(f[0], {0, 0});

        for (int j = 1; j <= m; j++) {
            update(f[j], f[j - 1] + Value{1, 0});
        }

        if (f[m] < opt) {
            opt = f[m];
            opt_i = i;
        }
    }

    Alignment result;
    result.range1 = {opt_i - opt.l + 1, opt_i + 1};
    result.range2 = {1, m + 1};
    result.loss = opt.t;

    return result;
}

auto Index::align(const BioSeq &s) -> Alignment {
    int n = s.size();
    std::priority_queue<State, std::vector<State>, Heuristic> q(Heuristic{n});
    tsl::robin_map<Key, int> best;

    auto probe = [&best, &q](const State &v) {
        auto it = best.find(v.key);

        if (it == best.end()) {
            best.insert({v.key, v.t});
            q.push(v);
        } else if (it->second > v.t) {
            it.value() = v.t;
            q.push(v);
        }
    };

    probe(State());

    State opt;
    size_t max_queue_size = 0;

    do {
        max_queue_size = std::max(max_queue_size, q.size());
        auto u = q.top();
        q.pop();

        if (u.t > best[u.key])
            continue;

        if (u.key.y == n) {
            opt = u;
            break;
        }

        auto [x, y] = u.key;

        probe({{x, y + 1}, u.t + FULL_COST, u.l});

        for (int c = 0; c < ALPHABET_SIZE; c++) {
            int z = m[x].transition[c];
            if (!z)
                continue;

            int cost = (c == CMAP[s[y + 1]] ? 0 : MISS_COST) + CHAR_COST;
            probe({{z, y + 1}, u.t + cost, u.l + 1});
            probe({{z, y}, u.t + FULL_COST, u.l + 1});
        }
    } while (!q.empty());

    AlignmentDebugInfo debug;
    debug.n_state_visited = best.size();
    debug.max_queue_size = max_queue_size;

    return Alignment{
        {opt.key.x, opt.l},
        opt.t, debug
    };
}

}
