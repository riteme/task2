#include <queue>
#include <unordered_map>

#include "core.hpp"

#include "tsl/robin_map.h"


namespace {

constexpr int KMER = 20;
constexpr int STEP = 3;
constexpr int MIN_BUCKET_SIZE = 850;
constexpr int NUM_SEQ = 2;
constexpr int MIN_THRESHOLD = 10;

}

namespace core {

auto Index::fuzzy_locate(const BioSeq &seq) const -> Location {
    int n = seq.size();

    std::string rev_seq = watson_crick_complement(*seq.internal);
    BioSeq s[NUM_SEQ] = {seq, rev_seq};

    int bucket_size = std::max(MIN_BUCKET_SIZE, n / 2);
    // tsl::robin_map<int, int> bucket[2];
    std::unordered_map<int, int> bucket[2];

    auto put = [&bucket, bucket_size](int i, int j) {
        j /= bucket_size;
        bucket[i][j]++;
    };

    auto probe = [&bucket](int i, int j) -> int {
        auto it = bucket[i].find(j);
        if (it != bucket[i].end())
            return it->second;
        return 0;
    };

    for (int i = 0; i < NUM_SEQ; i++) {
        for (int l = 1; l + KMER - 1 <= n; l += STEP) {
            auto t = align(s[i].take(l, l + KMER)).token;

            for (int j : rpset(t)) {
                put(i, j - t.len / 2);
            }
        }
    }

    int threshold = std::numeric_limits<int>::max();
    int max_score = std::numeric_limits<int>::min(), best_i = 0, best_j = 0;
    for (int i = 0; i < 2; i++) {
        for (auto &p : bucket[i]) {
            int j = p.first;
            int self = p.second;
            int prev = probe(i, j - 1);
            int succ = probe(i, j + 1);

            if (self * 2 < prev + succ)
                continue;

            int score = prev + self + succ;

            if (score > max_score) {
                int hi = std::max(self, std::max(prev, succ));
                int lo = std::min(self, std::min(prev, succ));
                int gap = hi - lo;
                threshold = lo - gap / 2;

                max_score = score;
                best_i = i;
                best_j = j;
            }
        }
    }

    threshold = std::max(threshold, MIN_THRESHOLD);

    int left = best_j - 1;
    while (probe(best_i, left - 1) >= threshold) {
        left--;
    }

    int right = best_j + 1;
    while (probe(best_i, right + 1) >= threshold) {
        right++;
    }

    // printf("bucket[%d][%d..%d] = {", best_i, left - 5, right + 5);
    // for (int k = left - 5; k <= right + 5; k++) {
    //     if (k == left)
    //         printf("[");
    //     printf("%d", bucket[best_i][k]);
    //     if (k == right)
    //         printf("]");
    //     if (k != right + 5)
    //         printf(" ");
    // }
    // printf("}\n");

    Location result;
    result.reversed = best_i == 0 ? false : true;
    result.left = std::max(1, left * bucket_size);
    result.right = std::min(size(), (right + 2) * bucket_size - 1);

    return result;
}

}
