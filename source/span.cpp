#include "index.hpp"
#include "numeric.hpp"


namespace {

constexpr int INF = 0x3f3f3f3f;

struct Record {
    int t, l1, l2;

    static auto zero() -> Record {
        return {0, 0, 0};
    }

    static auto max() -> Record {
        return {INF, 0, 0};
    }

    auto rate() const -> double {
        return 1.0 - double(t) / l2;
    }

    auto total_length() const -> int {
        return l1 + l2;
    }

    auto length_gap() const -> int {
        return std::abs(l1 - l2);
    }

    bool better_than(const Record &rhs) const {
        if (total_length() == rhs.total_length())
            return length_gap() < rhs.length_gap();
        return total_length() > rhs.total_length();
    }

    bool operator<(const Record &rhs) const {
        return t == rhs.t ? total_length() > rhs.total_length() : t < rhs.t;
    }

    auto operator+(const Record &rhs) const -> Record {
        return {t + rhs.t, l1 + rhs.l1, l2 + rhs.l2};
    }
};

template <typename T>
void update(T &dest, const T &value) {
    if (value < dest)
        dest = value;
}

}

namespace core {

template <typename TCompare, typename TOutput, bool Debug = false>
static inline auto _partial_span_impl(
    const BioSeq &s1, const BioSeq &s2,
    const TCompare &compare,
    const TOutput &output
) -> Alignment {
    constexpr int PENALTY = 3;

    constexpr auto MIN_DEVIATION = 1.0;
    constexpr auto DEVIATION_SCALE = 1.5;
    constexpr auto MIN_SLOPE = 0.8;

    constexpr int PUSH_STARTUP = 100;
    constexpr int PUSH_LOOKAFTER = 50;
    constexpr int SCAN_LOOKAHEAD = 30;
    constexpr int SCAN_LOOKAFTER = 50;
    constexpr int SCAN_DEVIATION_THRESHOLD = 20;
    constexpr int SCAN_MIN_SIZE = 32;

    int n = s1.size(), m = s2.size();

    std::vector<Record> f[2];
    for (int i = 0; i < 2; i++) {
        f[i].resize(m + 1);
        for (int j = 0; j <= m; j++) {
            f[i][j] = {j, 0, j};
        }
    }

    std::vector<Record> opt;
    opt.resize(m + 1, Record::max());

    for (int i = 1; i <= n; i++) {
        for (int j = m; j > 0; j--) {
            f[1][j] = f[1][j] + Record{1, 1, 0};
            update(f[1][j], f[0][j] + Record{1 + PENALTY, 1, 0});

            f[0][j] = Record::max();
            if (compare(i, j)) {
                update(f[0][j], f[0][j - 1] + Record{0, 1, 1});
                update(f[0][j], f[1][j - 1] + Record{0, 1, 1});
            }
        }

        f[0][0] = Record::max();
        f[1][0] = f[1][0] + Record{1, 1, 0};
        update(f[1][0], f[0][0] + Record{1 + PENALTY, 1, 0});

        for (int j = 1; j <= m; j++) {
            update(f[1][j], f[1][j - 1] + Record{1, 0, 1});
            update(f[1][j], f[0][j - 1] + Record{1 + PENALTY, 0, 1});
        }

        for (int j = 0; j <= m; j++) {
            update(opt[j], f[0][j]);
            update(opt[j], f[1][j]);
        }
    }

    std::vector<Vec2d> vs;
    vs.resize(m + 1);
    for (int j = 0; j <= m; j++) {
        vs[j] = Vec2d(j, opt[j].l1);
    }

    if (Debug) {
        printf("Show[ListPlot[{");
        for (int j = 1; j <= m; j++) {
            printf("{%d,%d}", j, opt[j].l1);
            if (j < m)
                printf(",");
        }
        puts("}]");
    }

    auto deviation = [](const Vec2d &line, const Vec2d &point) {
        auto y0 = line.k() * point.x + line.b();
        return std::abs(point.y - y0);
    };

    int r = PUSH_STARTUP, last_r = 0;
    Vec2d primary_line;
    while (last_r != r) {
        last_r = r;
        primary_line = linear_least_square(vs.begin(), vs.begin() + r, 8);

        auto max_dev = MIN_DEVIATION;
        for (int i = 0; 2 * i < r; i++) {
            max_dev = std::max(max_dev, deviation(primary_line, vs[i]));
        }

        if (Debug)
            fprintf(stderr, "max_dev=%.4lf\n", max_dev);
        max_dev *= DEVIATION_SCALE;

        for (int i = std::min(m, r + PUSH_LOOKAFTER); i >= r; ) {
            if (deviation(primary_line, vs[i]) < max_dev) {
                r = i + 1;
                i = std::min(m, r + PUSH_LOOKAFTER);
            } else
                i--;
        }

        if (Debug) {
            printf(
                ",Plot[Style[(%.16lf)x+(%.16lf),RGBColor[0,0,0,%.2lf]],{x,0,%d}]\n",
                primary_line.k(), primary_line.b(), (r == last_r ? 1.0 : 0.2), m
            );
            fprintf(stderr, "r=%d\n", r);
        }
    }

    std::vector<Vec2d> neighbors;

    if (primary_line.k() > MIN_SLOPE) {
        neighbors.reserve(SCAN_LOOKAFTER);
        for (int i = std::max(0, r - SCAN_LOOKAHEAD); i <= m && i <= r + SCAN_LOOKAFTER; i++) {
            if (std::abs(opt[r - 1].l1 - opt[i].l1) <= SCAN_DEVIATION_THRESHOLD)
                neighbors.push_back(vs[i]);
        }
    }

    if (Debug)
        fprintf(stderr, "neighbors.size()=%zu\n", neighbors.size());

    int corner = 0;
    if (neighbors.size() > SCAN_MIN_SIZE) {
        auto secondary_line = linear_least_square(neighbors.begin(), neighbors.end());
        auto [x, y] = line_intersection(primary_line, secondary_line);

        if (Debug) {
            printf(
                ",Plot[Style[(%.16lf)x+(%.16lf),Black],{x,0,%d}]\n",
                secondary_line.k(), secondary_line.b(), m
            );
            printf(
                ",Epilog->{Directive[RGBColor[1,0,0,0.5]],Line[{{0,%.16lf},{%d,%.16lf}}],Line[{{%.16lf,0},{%.16lf,%d}}]},ImageSize->Full]\n",
                y, m, y, x, x, n
            );
        }

        corner = std::max(0, std::min(m, static_cast<int>(std::round(x))));
    } else {
        if (Debug)
            puts(",ImageSize->Full]");
    }

    auto best = opt[corner];
    return output(best, best.l1, best.l2);
}

auto prefix_span(const BioSeq &s1, const BioSeq &s2) -> Alignment {
    return _partial_span_impl(
        s1, s2,
        [&](int i, int j) {
            return s1[i] == s2[j];
        },
        [&](const Record &opt, int i, int j) {
            Alignment result;
            result.range1 = {1, i + 1};
            result.range2 = {1, j + 1};
            result.loss = opt.t;
            return result;
        }
    );
}

auto suffix_span(const BioSeq &s1, const BioSeq &s2) -> Alignment {
    int n = s1.size(), m = s2.size();

    return _partial_span_impl(
        s1, s2,
        [&](int i, int j) {
            return s1[n - i + 1] == s2[m - j + 1];
        },
        [&](const Record &opt, int i, int j) {
            Alignment result;
            result.range1 = {n - i + 1, n + 1};
            result.range2 = {m - j + 1, m + 1};
            result.loss = opt.t;
            return result;
        }
    );
}

}
