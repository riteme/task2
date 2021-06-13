#pragma once

#include <cmath>

#include <vector>


namespace core {

template <typename T>
struct Vec2 {
    Vec2() : x(0.0), y(0.0) {}
    Vec2(T _x, T _y) : x(_x), y(_y) {}

    template <typename U>
    Vec2(const Vec2<U> &rhs) : x(rhs.x), y(rhs.y) {}

    template <typename U>
    auto operator=(const Vec2<U> &rhs) -> Vec2 & {
        x = rhs.x;
        y = rhs.y;
        return *this;
    }

    T x, y;

    auto len2() const -> T {
        return x * x + y * y;
    }

    auto len() const -> double {
        return hypot(x, y);
    }

    auto k() const -> double {
        return x;
    }

    auto k() -> double & {
        return x;
    }

    auto b() const -> double {
        return y;
    }

    auto b() -> double & {
        return y;
    }

    auto unit() const -> Vec2<double> {
        return Vec2<double>(*this) / len();
    }

    auto operator+(const Vec2 &rhs) const -> Vec2 {
        return {x + rhs.x, y + rhs.y};
    }

    auto operator-(const Vec2 &rhs) const -> Vec2 {
        return {x - rhs.x, y - rhs.y};
    }

    auto operator*(const T &rhs) const -> Vec2 {
        return {x * rhs, y * rhs};
    }

    auto operator/(const T &rhs) const -> Vec2 {
        return {x / rhs, y / rhs};
    }

    auto operator*(const Vec2 &rhs) const -> T {
        return x * rhs.x + y * rhs.y;
    }

    auto operator%(const Vec2 &rhs) const -> T {
        return x * rhs.y - y * rhs.x;
    }
};

using Vec2i = Vec2<int>;
using Vec2f = Vec2<float>;
using Vec2d = Vec2<double>;

template <typename TIterator>
concept Vec2dIterator = requires (TIterator it) {
    static_cast<Vec2d>(*it);
    static_cast<double>(it->x);
    static_cast<double>(it->y);
};

template <Vec2dIterator TIterator>
auto linear_least_square(TIterator beg, TIterator end, int n_reduce = 0) -> Vec2d {
    constexpr int N_THRESHOLD = 30;

    int n = 0;

    auto sx = 0.0, sy = 0.0, sxy = 0.0, sx2 = 0.0;
    for (auto it = beg; it != end; it++, n++) {
        sx += it->x;
        sy += it->y;
        sxy += it->x * it->y;
        sx2 += it->x * it->x;
    }

    auto k = (n * sxy - sx * sy) / (n * sx2 - sx * sx);

    if (n_reduce > 0) {
        auto b = (sy - k * sx) / n;

        auto dev = [k, b](auto it) {
            return std::abs(it->y - (k * it->x + b));
        };

        auto sdev = 0.0;
        for (auto it = beg; it != end; it++) {
            sdev += dev(it);
        }

        auto threshold = 2 * sdev / n;
        std::vector<Vec2d> vs;
        vs.reserve(n);
        for (auto it = beg; it != end; it++) {
            if (dev(it) <= threshold)
                vs.push_back(*it);
        }

        if (N_THRESHOLD <= vs.size() && vs.size() < n)
            return linear_least_square(vs.begin(), vs.end(), n_reduce - 1);
    }

    auto b = (sy - k * sx) / n;
    return {k, b};
}

auto line_intersection(const Vec2d &l1, const Vec2d &l2) -> Vec2d;

struct Decomposition {
    std::vector<int> start;
    double area;
};

auto french_stick_decompose(const std::vector<Vec2d> &vs, int K) -> Decomposition;

}
