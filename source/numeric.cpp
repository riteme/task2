#include <cstdio>

#include "numeric.hpp"


namespace core {

auto line_intersection(const Vec2d &l1, const Vec2d &l2) -> Vec2d {
    auto x = (l2.b() - l1.b()) / (l1.k() - l2.k());
    auto y = l1.k() * x + l1.b();
    return {x, y};
}

}
