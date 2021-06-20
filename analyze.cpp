#include <vector>
#include <sstream>
#include <fstream>
#include <functional>
#include <unordered_set>
#include <unordered_map>

#include "CLI11.hpp"

#include "core.hpp"


namespace {

constexpr int SNAP_DISTANCE = 200;
constexpr int MAX_SV_LENGTH = 1100;
constexpr int MIN_SV_LENGTH = 50;
constexpr auto INV_MIN_SCORE = 0.65;
constexpr int MAX_CONJECTION_LENGTH = 150;
constexpr auto MIN_CONJECTION_MATCH_RATE = 0.6;
constexpr auto MAX_TRA_DISCREPANCY = 20.0;
constexpr int LOCATOR_LENGTH = 100;
constexpr int EXTRA_LOCATOR_LENGTH = 256;
constexpr int SCAN_LENGTH = 1800;
constexpr auto LOCATOR_MIN_MATCH_RATE = 0.75;

// link type
enum class LType : int {
    INV, DUP, DEL, INS, TRA
};

auto to_string(const LType &type) -> const char * {
    switch (type) {
        case LType::INV: return "INV";
        case LType::DUP: return "DUP";
        case LType::DEL: return "DEL";
        case LType::INS: return "INS";
        case LType::TRA: return "TRA";
        default: return "(unknown)";
    }
}

struct Link {
    LType type;
    struct Endpoint *ep;
};

using AdjList = std::vector<Link>;

// endpoint type.
enum class EType {
    LEFT, RIGHT
};

struct Endpoint {
    std::string name;
    int pos1, pos2, len;
    AdjList adj;
    bool marked = false;

    bool empty() const {
        return adj.empty();
    }

    bool snap_to(const Endpoint &ep, int max_dist = SNAP_DISTANCE) {
        return snap_to(ep.pos1, max_dist);
    }

    bool snap_to(int pos, int max_dist = SNAP_DISTANCE) {
        return std::abs(pos - pos1) <= max_dist;
    }

    bool operator<(const Endpoint &rhs) const {
        return pos1 < rhs.pos1;
    }
};

// endpoint list.
using EList = std::vector<Endpoint>;
using ERefList = std::vector<Endpoint *>;

inline auto dist(const Endpoint &lp, const Endpoint &rp) -> int {
    return std::abs(lp.pos1 - rp.pos1);
}

inline void link(const LType &type, Endpoint &u, Endpoint &v) {
    u.adj.push_back({type, &v});
    v.adj.push_back({type, &u});
}

struct Range {
    int left, right;
    double score;

    auto length() const -> int {
        return right - left + 1;
    }
};

// range list.
using RList = std::vector<Range>;

struct Key {
    std::string name;
    EType type;

    bool operator==(const Key &) const = default;
};

// endpoint map.
using EMap = std::unordered_map<Key, EList>;

// range map.
using RMap = std::unordered_map<std::string, RList>;

}

template <>
struct std::hash<EType> {
    auto operator()(const EType &u) const -> size_t {
        return u == EType::LEFT ? 0 : size_t(-1L);
    }
};

template <>
struct std::hash<Key> {
    auto operator()(const Key &u) const -> size_t {
        static hash<std::string> H1;
        static hash<EType> H2;
        return H1(u.name) ^ H2(u.type);
    }
};

namespace {

template <typename TStream, typename TData>
requires requires (TStream &s, TData &x) {
    s >> x;
}
void take(TStream &s, TData &dst) {
    s >> dst;
}

void load_locate_file(const std::string &path, core::Dict &runs) {
    std::fstream fp(path);
    while (fp) {
        std::string line;
        std::getline(fp, line);
        if (line.empty())
            continue;

        std::stringstream s(line);

        std::string name, _;
        bool reversed;

        take(s, name);
        take(s, _);  // ref_name
        take(s, _);  // left
        take(s, _);  // right
        take(s, _);  // loss
        take(s, reversed);

        if (reversed) {
            auto ptr = runs.find(name);
            ptr->sequence = core::watson_crick_complement(ptr->sequence);
        }
    }
}

auto load_dump_file(const std::string &path) {
    EMap emap;
    RMap rmap;

    std::fstream fp(path);
    while (fp) {
        std::string line;
        std::getline(fp, line);
        if (line.empty())
            continue;

        std::stringstream s(line);

        std::string run_name, ref_name;
        take(s, run_name);
        take(s, ref_name);

        Endpoint lp, rp;
        lp.name = run_name;
        rp.name = run_name;
        take(s, lp.pos1);
        take(s, lp.pos2);
        take(s, lp.len);
        take(s, rp.pos1);
        take(s, rp.pos2);
        take(s, rp.len);

        double inv_score;
        take(s, inv_score);

        if (lp.pos1 > 0)
            emap[{ref_name, EType::LEFT}].push_back(lp);
        if (rp.pos1 > 0)
            emap[{ref_name, EType::RIGHT}].push_back(rp);
        if (lp.pos1 > 0 && rp.pos1 > 0 && lp.pos1 < rp.pos1)
            rmap[ref_name].push_back({lp.pos1, rp.pos1, inv_score});
    }

    return std::make_tuple(emap, rmap);
}

}

int main(int argc, char *argv[]) {
    std::string ref_file, runs_file, dump_file, locate_file;

    CLI::App args;
    args.add_option("-r", ref_file)->required();
    args.add_option("-l", runs_file)->required();
    args.add_option("-p", locate_file)->required();
    args.add_option("-d", dump_file)->required();
    CLI11_PARSE(args, argc, argv);

    /**
     * load data.
     */

    core::Dict refs, runs;

    refs.load_file(ref_file);
    printf("loaded \"%s\".\n", ref_file.data());

    runs.load_file(runs_file);
    runs.build_index();
    printf("loaded \"%s\".\n", runs_file.data());

    load_locate_file(locate_file, runs);

    EMap emap;
    RMap rmap;
    std::tie(emap, rmap) = load_dump_file(dump_file);
    printf("loaded \"%s\".\n", dump_file.data());

    /**
     * probe special SVs.
     */

    auto probe_inv = [](EList &L, EList &R, RList &rs) {
        for (auto [l, r, score] : rs) {
            if (score < INV_MIN_SCORE)
                continue;

            for (auto &lp : L) for (auto &rp : R) {
                if (lp.snap_to(l) && rp.snap_to(r))
                    link(LType::INV, lp, rp);
            }
        }
    };

    auto probe_del_and_dup = [&](EList &L, EList &R) {
        for (auto &lp : L) for (auto &rp : R) {
            if (dist(lp, rp) < MIN_SV_LENGTH ||
                dist(lp, rp) > MAX_SV_LENGTH)
                continue;

            auto &seq1 = runs.find(lp.name)->sequence;
            auto &seq2 = runs.find(rp.name)->sequence;
            int size1 = seq1.size();
            int size2 = seq2.size();

            int left_len = std::min(
                MAX_CONJECTION_LENGTH,
                std::min(lp.pos2, rp.pos2)
            );
            int right_len = std::min(
                MAX_CONJECTION_LENGTH,
                std::min(size1 - lp.pos2, size2 - rp.pos2)
            );
            int len = left_len + right_len;

            auto slice1 = core::BioSeq(seq1,
                std::max(1, lp.pos2 - left_len + 1),
                std::min(size1 + 1, lp.pos2 + right_len)
            );
            auto slice2 = core::BioSeq(seq2,
                std::max(1, rp.pos2 - left_len + 1),
                std::min(size2 + 1, rp.pos2 + right_len)
            );

            int loss = core::full_align(slice1, slice2);

            auto rate = 1 - double(loss) / len;

            // if (rate > 0.5)
            //     printf("rate=%.4lf, lp=%d, rp=%d\n", rate, lp.pos1, rp.pos1);

            if (rate >= MIN_CONJECTION_MATCH_RATE) {
                if (lp < rp)
                    link(LType::DEL, lp, rp);
                else
                    link(LType::DUP, lp, rp);
            }
        }
    };

    auto probe_ins = [](EList &L, EList &R) {
        for (auto &lp : L) for (auto &rp : R) {
            if (lp.snap_to(rp, MIN_SV_LENGTH))
                link(LType::INS, lp, rp);
        }
    };

    for (auto &[name, rs] : rmap) {
        auto &L = emap[{name, EType::LEFT}];
        auto &R = emap[{name, EType::RIGHT}];

        probe_inv(L, R, rs);
        probe_del_and_dup(L, R);
        probe_ins(L, R);
    }

    /**
     * aggregate and output.
     */

    auto reset = [&]() {
        for (auto &[_, es] : emap) {
            for (auto &ep : es) {
                ep.marked = false;
            }
        }
    };

    auto collect = [](const LType &type, Endpoint *x) {
        std::function<void(Endpoint *, ERefList &, ERefList &)> dfs;
        dfs = [&dfs, type](Endpoint *u, ERefList &L, ERefList &R) {
            if (u->marked)
                return;
            u->marked = true;

            L.push_back(u);
            for (auto e : u->adj) {
                if (e.type == type)
                    dfs(e.ep, R, L);
            }
        };

        ERefList L, R;
        dfs(x, L, R);
        return std::make_tuple(L, R);
    };

    auto accumulate = [](double &sum, int &count, const ERefList &es) {
        for (auto &ep : es) {
            sum += ep->pos1;
            count++;
        }
    };

    auto dump_normal = [&](
        const char *op,
        const std::string &name,
        ERefList L, ERefList R
    ) {
        auto sum = 0.0;
        int count = 0;
        accumulate(sum, count, L);
        int left = std::round(sum / count);

        sum = 0.0;
        count = 0;
        accumulate(sum, count, R);
        int right = std::round(sum / count);

        if (right < left)
            std::swap(left, right);

        fprintf(stderr, "%s %s %d %d\n", op, name.data(), left, right);
    };

    auto dump_ins = [&](
        const char *op,
        const std::string &name,
        ERefList L, ERefList R
    ) {
        auto sum = 0.0;
        int count = 0;
        accumulate(sum, count, L);
        accumulate(sum, count, R);
        int left = std::round(sum / count);

        sum = 0.0;
        count = 0;
        std::unordered_set<std::string> mark;

        for (auto &ep : L) {
            if (ep->len > MAX_SV_LENGTH)
                continue;

            sum += ep->len;
            count++;
            mark.insert(ep->name);
        }

        for (auto &ep : R) {
            if (ep->len > MAX_SV_LENGTH)
                continue;

            int scale = mark.count(ep->name) ? 4 : 1;
            sum += ep->len * scale;
            count += scale;
        }

        int right = left + std::round(sum / count);

        fprintf(stderr, "%s %s %d %d\n", op, name.data(), left, right);
    };

    auto dump = [&]<typename TDumpFn>(const LType &type, const TDumpFn &dump_fn) {
        reset();

        for (auto &e : refs) {
            for (auto &ep : emap[{e.name, EType::LEFT}]) {
                if (!ep.marked) {
                    auto [L, R] = collect(type, &ep);
                    if (!L.empty() && !R.empty())
                        dump_fn(to_string(type), e.name, L, R);
                }
            }
        }
    };

    dump(LType::INV, dump_normal);
    dump(LType::DEL, dump_normal);
    dump(LType::DUP, dump_normal);
    dump(LType::INS, dump_ins);

    /**
     * naïve TRA pairing & dumping.
     */

    struct Position {
        double pos;
        bool marked = false;

        void mark() {
            marked = true;
        }

        auto to_int() const -> int {
            return static_cast<int>(std::round(pos));
        }
    };

    using PosList = std::vector<Position>;
    using PosMap = std::unordered_map<Key, PosList>;

    auto compact = [&] {
        PosMap pmap;

        for (auto &[key, es] : emap) {
            auto &compacted = pmap[key];

            std::vector<double> list;
            list.reserve(es.size());
            for (auto &ep : es) {
                list.push_back(ep.pos1);
            }

            std::sort(list.begin(), list.end());

            for (auto i = list.begin(); i != list.end(); ) {
                auto j = i, k = std::next(i);
                while (k != list.end() && std::abs(*j - *k) <= SNAP_DISTANCE) {
                    j = k;
                    k++;
                }

                auto sum = 0.0;
                for (auto p = i; p != k; p++) {
                    sum += *p;
                }

                compacted.push_back({sum / (k - i), false});

                i = k;
            }
        }

        return pmap;
    };

    auto dump_tra = [&](PosMap &pmap) {
        for (auto &e1 : refs)
        for (auto &l1 : pmap[{e1.name, EType::LEFT}])
        for (auto &r1 : pmap[{e1.name, EType::RIGHT}]) {
            auto len1 = r1.pos - l1.pos;
            if (len1 < MIN_SV_LENGTH || len1 > MAX_SV_LENGTH)
                continue;

            for (auto &e2 : refs) if (e2.name > e1.name)
            for (auto &l2 : pmap[{e2.name, EType::LEFT}])
            for (auto &r2 : pmap[{e2.name, EType::RIGHT}]) {
                auto len2 = r2.pos - l2.pos;
                if (len2 < MIN_SV_LENGTH || len2 > MAX_SV_LENGTH)
                    continue;

                if (std::abs(len1 - len2) <= MAX_TRA_DISCREPANCY) {
                    int left1 = l1.to_int(), right1 = r1.to_int();
                    int left2 = l2.to_int(), right2 = r2.to_int();

                    l1.marked = true;
                    l2.marked = true;
                    r1.marked = true;
                    r2.marked = true;

                    fprintf(stderr,
                        "TRA %s %d %d %s %d %d\n",
                        e1.name.data(), left1, right1,
                        e2.name.data(), left2, right2
                    );
                }
            }
        }
    };

    auto pmap = compact();
    dump_tra(pmap);

    /**
     * extra dumping.
     */

    auto dump_extra_del_and_dup = [&] {
        for (auto &e : refs) {
            for (auto &lp : emap[{e.name, EType::LEFT}]) {
                auto &run = *runs.find(lp.name);
                if (!lp.empty() || lp.pos2 >= run.sequence.size() - LOCATOR_LENGTH)
                    continue;

                int len = std::min(EXTRA_LOCATOR_LENGTH, int(run.sequence.size()) - lp.pos2);
                auto t = core::BioSeq(run.sequence, lp.pos2 + 1, lp.pos2 + len);

                // DEL
                int right = std::min(int(e.sequence.size()), lp.pos1 + SCAN_LENGTH);
                auto s = core::BioSeq(e.sequence, lp.pos1 + 1, right);

                auto result = core::local_align(s, t);
                int pos = lp.pos1 + result.range1.begin;

                if (result.match_rate2() > LOCATOR_MIN_MATCH_RATE &&
                    std::abs(pos - lp.pos1) > MIN_SV_LENGTH) {
                    fprintf(stderr,
                        "DEL %s %d %d\n",
                        e.name.data(), lp.pos1, pos
                    );
                }

                // DUP
                int left = std::max(1, lp.pos1 - SCAN_LENGTH);
                s = core::BioSeq(e.sequence, left, lp.pos1);

                result = core::local_align(s, t);
                pos = left + result.range1.begin;

                if (result.match_rate2() > LOCATOR_MIN_MATCH_RATE &&
                    std::abs(pos - lp.pos1) > MIN_SV_LENGTH) {
                    fprintf(stderr,
                        "DUP %s %d %d\n",
                        e.name.data(), pos, lp.pos1
                    );
                }
            }

            for (auto &rp : emap[{e.name, EType::RIGHT}]) {
                if (!rp.empty() || rp.pos2 <= LOCATOR_LENGTH)
                    continue;

                auto &run = *runs.find(rp.name);
                int len = std::min(EXTRA_LOCATOR_LENGTH, rp.pos2 - 1);
                auto t = core::BioSeq(run.sequence, rp.pos2 - len, rp.pos2);

                // DEL
                int left = std::max(1, rp.pos1 - SCAN_LENGTH);
                auto s = core::BioSeq(e.sequence, left, rp.pos1);

                auto result = core::local_align(s, t);
                int pos = left + result.range1.end;

                if (result.match_rate2() > LOCATOR_MIN_MATCH_RATE &&
                    std::abs(rp.pos1 - pos) > MIN_SV_LENGTH) {
                    fprintf(stderr,
                        "DEL %s %d %d\n",
                        e.name.data(), pos, rp.pos1
                    );
                }

                // DUP
                int right = std::min(int(e.sequence.size()), rp.pos1 + SCAN_LENGTH);
                s = core::BioSeq(e.sequence, rp.pos1 + 1, right);

                result = core::local_align(s, t);
                pos = rp.pos1 + result.range1.end;

                if (result.match_rate2() > LOCATOR_MIN_MATCH_RATE &&
                    std::abs(rp.pos1 - pos) > MIN_SV_LENGTH) {
                    fprintf(stderr,
                        "DUP %s %d %d\n",
                        e.name.data(), rp.pos1, pos
                    );
                }
            }
        }
    };

    auto dump_extra_inv = [&]() {
        for (auto &e : refs)
        for (auto &l : pmap[{e.name, EType::LEFT}])
        for (auto &r : pmap[{e.name, EType::RIGHT}]) {
            auto len = r.pos - l.pos;
            if (!l.marked && !r.marked &&
                len >= MIN_SV_LENGTH &&
                len <= MAX_SV_LENGTH) {
                int left = l.to_int();
                int right = r.to_int();

                fprintf(stderr,
                    "INV %s %d %d\n",
                    e.name.data(), left, right
                );
            }
        }
    };

    dump_extra_del_and_dup();
    dump_extra_inv();

    return 0;
}
