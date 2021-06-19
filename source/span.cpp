#include <tuple>

#include "index.hpp"
#include "numeric.hpp"


namespace {

using namespace core;

constexpr int INF = 0x3f3f3f3f;
constexpr auto MIN_SLOPE = 0.8;

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

template <Vec2dIterator TIterator>
auto range_slope(TIterator beg, const TIterator &end) -> double {
    int n = 0;
    auto min_x = std::numeric_limits<double>::max();
    auto max_x = std::numeric_limits<double>::min();
    auto min_y = min_x, max_y = max_x;

    for (auto it = beg; it != end; it++, n++) {
        min_x = std::min(min_x, it->x);
        max_x = std::max(max_x, it->x);
        min_y = std::min(min_y, it->y);
        max_y = std::max(max_y, it->y);
    }

    if (n < 2)
        return 1.0;
    else if (n < 5)
        return (max_y - min_y) / 10;
    else
        return (max_y - min_y) / std::max(0.1, max_x - min_x);
}

template <typename T, typename THeightFn>
requires requires (T x, THeightFn fn) {
    static_cast<int>(fn(x));
}
auto trim_outliers(const std::vector<T> &vs, const THeightFn &height) -> std::vector<bool> {
    constexpr int INNER_THRESHOLD = 50;

    int n = vs.size();

    // x: index; y: corresponding height
    std::vector<Vec2i> f, bucket;
    f.resize(n);
    bucket.push_back({-1, std::numeric_limits<int>::min()});

    for (int i = 0; i < n; i++) {
        int y = height(vs[i]);
        int j = std::upper_bound(
            bucket.begin(), bucket.end(), Vec2i{-1, y},
            [](const Vec2i &u, const Vec2i &v) {
                return u.y < v.y;
            }
        ) - bucket.begin();

        f[i] = {bucket[j - 1].x, j};

        if (j == bucket.size())
            bucket.push_back({i, y});
        else if (bucket[j].y >= y)
            bucket[j] = {i, y};
    }

    std::vector<bool> mark;
    mark.resize(n);

    for (int i = bucket.back().x; i != -1; i = f[i].x) {
        mark[i] = true;
    }

    std::vector<bool> new_mark;
    new_mark.resize(n);

    int last = -1;
    for (int i = 0; i < n; i++) {
        if (mark[i])
            last = i;
        else if (last != -1 && std::abs(height(vs[last]) - height(vs[i])) <= INNER_THRESHOLD)
            new_mark[i] = true;
    }

    last = -1;
    for (int i = n - 1; i >= 0; i--) {
        if (mark[i])
            last = i;
        else if (last != -1 && std::abs(height(vs[last]) - height(vs[i])) <= INNER_THRESHOLD)
            new_mark[i] = true;
    }

    for (int i = 0; i < n; i++) {
        if (new_mark[i])
            mark[i] = true;
    }

    return mark;
}

auto decompose(std::vector<Vec2d> vs) -> Decomposition {
    constexpr int MIN_SLICE_LEN = 45;
    constexpr auto MAX_SLOPE = 9.5;
    constexpr auto SLOPE_DEVIATION_THRESHOLD = 0.1;
    constexpr int TAIL_CUT_MAX_LENGTH = 25;

    int K = 3;

    int offset = 0;
    Decomposition result;
    while (K > 0) {
        auto last_size = vs.size();
        result = french_stick_decompose(vs, K);

        // // experimental: head trim
        // if (result.slices.size() > 1) {
        //     auto &slice = result.slices[0];
        //     auto slope = range_slope(vs.begin(), vs.begin() + slice.end);
        //     fprintf(stderr, "slope=%.4lf\n", slope);
        //     if (slope < MIN_SLOPE) {
        //         offset = slice.end;
        //         vs.erase(vs.begin(), vs.begin() + offset);
        //         continue;
        //     }
        // }

        int fail_count = 0;
        for (auto &slice : result.slices) {
            if (slice.length() >= MIN_SLICE_LEN)
                continue;

            fail_count++;

            // tail erase
            bool do_erase = slice.end >= vs.size();

            // short cut
            if (!do_erase && slice.length() > 1) {
                auto slope = range_slope(vs.begin() + slice.begin, vs.begin() + slice.end);
                // fprintf(stderr, "slope=%.4lf\n", slope);
                do_erase = slope > MAX_SLOPE;
            }

            if (do_erase) {
                vs.erase(vs.begin() + slice.begin, vs.end());
                break;
            }
        }

        // slice meld
        if (fail_count == 0 && result.slices.size() > 1) {
            auto &s1 = result.slices[0];
            auto &s2 = result.slices[1];

            auto k1 = linear_least_square(vs.begin() + s1.begin, vs.begin() + s1.end).k();
            auto k2 = linear_least_square(vs.begin() + s2.begin, vs.begin() + s2.end).k();
            if (std::abs(k1 - k2) <= SLOPE_DEVIATION_THRESHOLD)
                fail_count++;
        }

        // tail cut
        if (fail_count == 0 && result.slices.size() > 1) {
            auto &s = result.slices[1];
            int len = std::min(s.length() / 2, TAIL_CUT_MAX_LENGTH);
            auto left = vs.begin() + s.end - len;
            auto slope = range_slope(left, vs.begin() + s.end);
            // fprintf(stderr, "slope=%.4lf\n", slope);
            if (slope > MAX_SLOPE) {
                vs.erase(left, vs.end());
                if (K == 3)
                    fail_count++;
            }
        }

        if (fail_count == 0 && vs.size() == last_size)
            break;

        K -= fail_count;
    }

    for (auto &slice : result.slices) {
        slice.begin += offset;
        slice.end += offset;
    }

    return result;
}

}

namespace core {

template <
    typename TFactory,
    bool Debug = false
>
static inline auto _partial_span_impl(
    const BioSeq &s1, const BioSeq &s2,
    const TFactory &factory,
    int offset = 0,
    bool enable_correlation = true
) -> Alignment {
    constexpr int PENALTY = 3;
    constexpr int N_REDUCE = 8;
    constexpr int OFFSET_THRESHOLD = 10;
    constexpr int LOCATOR_LENGTH = 100;

    auto output = factory(offset);

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
            if (s1[offset + i] == s2[j]) {
                update(f[0][j], f[0][j - 1] + Record{0, 1, 1});
                update(f[0][j], f[1][j - 1] + Record{0, 1, 1});
            }
        }

        f[1][0] = f[1][0] + Record{1, 1, 0};
        update(f[1][0], f[0][0] + Record{1 + PENALTY, 1, 0});
        f[0][0] = Record::max();

        for (int j = 1; j <= m; j++) {
            update(f[1][j], f[1][j - 1] + Record{1, 0, 1});
            update(f[1][j], f[0][j - 1] + Record{1 + PENALTY, 0, 1});
        }

        for (int j = 0; j <= m; j++) {
            update(opt[j], f[0][j]);
            update(opt[j], f[1][j]);
        }
    }

    auto mark = trim_outliers(opt, [](const Record &u) {
        return u.l1;
    });

    std::vector<Vec2d> vs;
    vs.reserve(m + 1);
    for (int j = 0; j <= m; j++) {
        if (mark[j])
            vs.push_back(Vec2d(opt[j].l2, opt[j].l1));
    }

    auto decomp = decompose(vs);
    int corner = int(vs[decomp.slices[0].end - 1].x + 0.5);
    if (corner > m)
        corner = m;
    if (corner < 0)
        corner = 0;

    if (Debug) {
        fprintf(stderr, "decomp.area = %.4lf\n", decomp.area);
        fprintf(stderr, "decomp.slices = ");
        for (auto &slice : decomp.slices) {
            fprintf(stderr, "[%d, %d) ", slice.begin, slice.end);
        }
        fprintf(stderr, "\n");

        const char *colors[] = {"RGBColor[1,0,0,0.1]", "RGBColor[0,1,0,0.1]", "RGBColor[0,0,1,0.1]"};
        printf("ListPlot[{");
        for (int i = 0; i < decomp.slices.size(); i++) {
            auto &slice = decomp.slices[i];

            printf("Style[{");
            for (int j = slice.begin; j < slice.end; j++) {
                printf("{%.1lf,%.1lf}", vs[j].x, vs[j].y);
                if (j + 1 < slice.end)
                    printf(",");
            }
            printf("},%s]\n", colors[i]);
            if (i + 1 < decomp.slices.size())
                printf(",");
        }

        int last = decomp.slices.back().end;
        if (last < vs.size()) {
            printf(",Style[{");
            for (int i = last; i < vs.size(); i++) {
                printf("{%.1lf,%.1lf}", vs[i].x, vs[i].y);
                if (i + 1 < vs.size())
                    printf(",");
            }
            printf("},RGBColor[0,0,0,0.1]]\n");
        }

        printf(",Style[{{0,0}");
        for (int j = 0; j <= m; j++) {
            if (!mark[j])
                printf(",{%d,%d}", opt[j].l2, opt[j].l1);
        }
        printf("},Orange]\n");

        int x = decomp.slices.size() > 1 ? opt[corner].l2 : m;
        int y = opt[x].l1;

        printf(
            "},Epilog->{Directive[Black],Line[{{0,%d},{%d,%d}}],Line[{{%d,0},{%d,%d}}]}",
            y, x, y, x, x, y
        );
        printf(",ImageSize->Full,PlotRange->All]\n");
    }

    auto slope = linear_least_square(vs.begin(), vs.begin() + decomp.slices[0].end, N_REDUCE).k();
    // fprintf(stderr, "slope=%.4lf\n", slope);
    bool slope_notify = slope < MIN_SLOPE;

    if (decomp.slices.size() <= 1 || slope_notify) {
        if (slope_notify && enable_correlation) {
            auto alignment = local_align(
                s1, s2.take(1, std::min(s2.size(), LOCATOR_LENGTH) + 1)
            );

            int offset = alignment.range1.begin;
            printf("warn: triggered correlation: offset=%d\n", offset);

            if (offset > OFFSET_THRESHOLD)
                return _partial_span_impl(s1, s2, factory, offset, false);
        }

        corner = 0;
    }

    auto best = opt[corner];
    auto result = output(best, best.l1, best.l2);
    result.mark = slope_notify;
    return result;
}

auto prefix_span(const BioSeq &s1, const BioSeq &s2) -> Alignment {
    return _partial_span_impl(
        s1, s2,
        [&](int offset) {
            return [&, offset](const Record &opt, int i, int j) {
                Alignment result;
                result.range1 = {1, offset + i + 1};
                result.range2 = {1, j + 1};
                result.loss = opt.t;
                return result;
            };
        }
    );
}

auto suffix_span(const BioSeq &_s1, const BioSeq &_s2) -> Alignment {
    int n = _s1.size(), m = _s2.size();

    auto c1 = _s1.clone(), c2 = _s2.clone();
    std::reverse(c1.begin(), c1.end());
    std::reverse(c2.begin(), c2.end());

    auto s1 = BioSeq(c1), s2 = BioSeq(c2);

    return _partial_span_impl(
        s1, s2,
        [&](int offset) {
            return [&, offset](const Record &opt, int i, int j) {
                Alignment result;
                result.range1 = {n - offset - i + 1, n + 1};
                result.range2 = {m - j + 1, m + 1};
                result.loss = opt.t;
                return result;
            };
        }
    );
}

}
