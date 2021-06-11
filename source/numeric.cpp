#include <cstdio>

#include "numeric.hpp"


namespace core {

auto bend_detect(const std::vector<Vec2d> &vs) -> int {
    constexpr int SPAN = 30;
    constexpr int LEAST_SQUARE_REDUCE = 8;
    constexpr double MAX_SLOPE = 2.0;
    constexpr double MIN_SLOPE = -0.2;

    std::vector<double> ks;
    ks.resize(vs.size());
    for (int i = SPAN; i + SPAN < vs.size(); i++) {
        ks[i] = linear_least_square(
            vs.begin() + i - SPAN,
            vs.begin() + i + SPAN + 1,
            LEAST_SQUARE_REDUCE
        );

        if (ks[i] > MAX_SLOPE)
            ks[i] = MAX_SLOPE;
        else if (ks[i] < MIN_SLOPE)
            ks[i] = MIN_SLOPE;
    }

    printf("ListPlot[{");
    for (int i = SPAN; i + SPAN < vs.size(); i++) {
        printf("{%d,%.4lf}", i,ks[i]);
        if (i + SPAN + 1 < vs.size())
            printf(",");
    }
    puts("},PlotRange -> All]");

    for (int i = SPAN; i + SPAN < vs.size(); i++) {
        if (ks[i] > 0.15)
            continue;

        int j = i;
        while (j >= SPAN && ks[j] < 0.7) {
            j--;
        }

        int k = i;
        while (k + SPAN < vs.size()) {
            bool ok = false;
            for (int p = k; p + SPAN < vs.size() && p < k + 32; p++) {
                ok |= ks[p] < 0.7 || ks[p] > 1.5;
            }

            if (ok)
                k++;
            else
                break;
        }

        printf("? j=%d, i=%d, k=%d, k-j=%d\n", j, i, k, k-j);
        if (k - j - 1 > 60)
            return (i + j) / 2;
    }

    return 0;
}

}
