#include <cstdio>

#include <queue>

#include "tsl/robin_map.h"

#include "index.hpp"


namespace {

constexpr int INF = 0x3f3f3f3f;

constexpr int MISS_COST = 10;
constexpr int CHAR_COST = 1;
constexpr int FULL_COST = MISS_COST + CHAR_COST;
constexpr int H_VALUE = 5;

struct Key {
    int x = 1, y = 0;

    bool operator==(const Key &rhs) const = default;

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

template <typename T>
void update(T &dest, const T &value) {
    if (value < dest)
        dest = value;
}

}

template <>
struct std::hash<Key> {
    auto operator()(const Key &z) const -> size_t {
        return ((z.x << 15) ^ z.y) * 0x19260817;
    }
};

namespace core {

auto full_align(const BioSeq &s1, const BioSeq &s2) -> int {
    int n = s1.size(), m = s2.size();
    std::vector<int> f;

    f.resize(m + 1);
    for (int j = 0; j <= m; j++) {
        f[j] = j;
    }

    for (int i = 1; i <= n; i++) {
        for (int j = m; j > 0; j--) {
            f[j] += 1;
            if (s1[i] == s2[j])
                f[j] = std::min(f[j], f[j - 1]);
        }

        f[0] += 1;

        for (int j = 1; j <= m; j++) {
            f[j] = std::min(f[j], f[j - 1] + 1);
        }
    }

    return f[m];
}

auto local_align(const BioSeq &s1, const BioSeq &s2) -> Alignment {
    constexpr int PENALTY = 3;

    struct Value {
        // d: (length in s1) - (length in s2) => n - m
        int t, d, loss;

        static auto zero() -> Value {
            return {0, 0, 0};
        }

        static auto max() -> Value {
            return {INF, INF, INF};
        }

        bool operator<(const Value &rhs) const {
            return t == rhs.t ? std::abs(d) < std::abs(rhs.d) : t < rhs.t;
        }

        auto operator+(const Value &rhs) const -> Value {
            return {t + rhs.t, d + rhs.d, loss + rhs.loss};
        }
    };

    int n = s1.size(), m = s2.size();

    std::vector<Value> f[2];
    for (int i = 0; i < 2; i++) {
        f[i].resize(m + 1);
        for (int j = 0; j <= m; j++) {
            f[i][j] = {j, -j, j};
        }
    }

    auto opt = Value::max();
    int opt_i = 0;
    for (int i = 1; i <= n; i++) {
        for (int j = m; j > 0; j--) {
            f[1][j] = f[1][j] + Value{1, +1, 1};
            update(f[1][j], f[0][j] + Value{1 + PENALTY, +1, 1});

            f[0][j] = Value::max();
            if (s1[i] == s2[j]) {
                update(f[0][j], f[0][j - 1] + Value{0, +0, 0});
                update(f[0][j], f[1][j - 1] + Value{0, +0, 0});
            }
        }

        f[1][0] = f[1][0] + Value{1, +1, 1};
        update(f[1][0], f[0][0] + Value{1 + PENALTY, +1, 1});
        f[0][0] = {0, 0, 0};

        for (int j = 1; j <= m; j++) {
            update(f[1][j], f[1][j - 1] + Value{1, -1, 1});
            update(f[1][j], f[0][j - 1] + Value{1 + PENALTY, -1, 1});
        }

        if (f[0][m] < opt) {
            opt = f[0][m];
            opt_i = i;
        }
        if (f[1][m] < opt) {
            opt = f[1][m];
            opt_i = i;
        }
    }

    Alignment result;
    int len = m + opt.d;
    result.range1 = {opt_i - len + 1, opt_i + 1};
    result.range2 = {1, m + 1};
    result.loss = opt.loss;

    return result;
}

auto Index::align(const BioSeq &s) const -> Alignment {
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
