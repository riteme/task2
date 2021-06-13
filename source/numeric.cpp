#include <cstdio>
#include <cassert>

#include <algorithm>
#include <functional>

#include "numeric.hpp"


namespace {

using namespace core;

auto last_edge(std::vector<Vec2d> &vs) -> double {
    if (vs.size() < 2)
        return 0;

    auto &p = vs[vs.size() - 1];
    auto &q = vs[vs.size() - 2];
    return (p % q) / 2;
};

constexpr bool Upper = true;
constexpr bool Lower = false;

template <bool Upper>
auto push_into(std::vector<Vec2d> &vs, const Vec2d &p) -> double {
    auto sum = 0.0;
    while (vs.size() > 1) {
        auto &q = vs[vs.size() - 1];
        auto &r = vs[vs.size() - 2];
        auto c = (q - p) % (r - p);

        bool pred = Upper ? c <= 0 : c >= 0;
        if (pred) {
            sum += last_edge(vs);
            vs.pop_back();
        } else
            break;
    }

    vs.push_back(p);
    return sum;
};

template <typename U>
concept DoubleIterator = requires (U u) {
    *u = 0.0;
};

template <typename T, typename U>
requires Vec2dIterator<T> && DoubleIterator<U>
void progressive_convex_hull(T beg, const T &end, U dest) {
    auto sum = 0.0;
    std::vector<Vec2d> upper, lower;

    for (auto it = beg; it != end; it++, dest++) {
        sum += push_into<Upper>(upper, *it);
        sum -= push_into<Lower>(lower, *it);
        sum += last_edge(lower);
        sum -= last_edge(upper);
        *dest = std::abs(sum);
    }
}

void vector_sqrt(std::vector<double> &vs) {
    constexpr auto BEND_COEFFICIENT = 0.45;

    for (auto &v : vs) {
        // v = std::sqrt(v);
        v = std::pow(v, BEND_COEFFICIENT);
    }
}

}

namespace core {

auto line_intersection(const Vec2d &l1, const Vec2d &l2) -> Vec2d {
    auto x = (l2.b() - l1.b()) / (l1.k() - l2.k());
    auto y = l1.k() * x + l1.b();
    return {x, y};
}

auto french_stick_decompose(const std::vector<Vec2d> &vs, int K) -> Decomposition {
    assert(K > 0);
    assert(std::is_sorted(vs.begin(), vs.end(), [](const Vec2d &u, const Vec2d &v) {
        return u.x < v.x;
    }));

    int n = vs.size();

    std::vector<double> suffix;
    suffix.resize(n);
    progressive_convex_hull(vs.rbegin(), vs.rend(), suffix.rbegin());
    vector_sqrt(suffix);

    std::function<Decomposition(int, int)> _decompose_impl;
    _decompose_impl = [n, &vs, &suffix, &_decompose_impl]
    (int K, int beg) -> Decomposition {
        if (K == 1)
            return {{{beg, n}}, suffix[beg]};

        int m = n - beg;

        std::vector<double> prefix;
        prefix.resize(m);
        progressive_convex_hull(vs.begin() + beg, vs.end(), prefix.begin());
        vector_sqrt(prefix);

        auto opt = Decomposition{{}, std::numeric_limits<double>::max()};
        int opt_i = 0;
        for (int i = 0; i + K <= m; i++) {
            if (prefix[i] > opt.area)
                break;

            auto subopt = _decompose_impl(K - 1, beg + i + 1);
            auto new_area = prefix[i] + subopt.area;
            if (opt.area > new_area) {
                opt_i = i;
                opt.area = new_area;
                opt.slices = std::move(subopt.slices);
            }
        }

        opt.slices.push_back({beg, beg + opt_i + 1});
        return opt;
    };

    auto result = _decompose_impl(K, 0);
    std::reverse(result.slices.begin(), result.slices.end());
    return result;
}

}
