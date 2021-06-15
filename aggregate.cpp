#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "CLI11.hpp"

#include "core.hpp"


enum class EndpointType {
    LEFT, RIGHT
};

struct Key {
    std::string name;
    EndpointType type;

    bool operator==(const Key &) const = default;
};

template <>
struct std::hash<Key> {
    auto operator()(const Key &u) const -> size_t {
        static hash<std::string> S;
        auto H = S(u.name);
        return (u.type == EndpointType::LEFT ? 0 : 0xffffffff) ^ H;
    }
};

template <typename T>
auto get(std::stringstream &stream) -> T {
    T data;
    stream >> data;
    return data;
}

struct UnionFind {
    std::vector<int> parent;

    void reset(int n) {
        parent.resize(n);
        for (int i = 0; i < n; i++) {
            parent[i] = i;
        }
    }

    auto root(int x) -> int {
        return parent[x] == x ? x : parent[x] = root(parent[x]);
    }

    void link(int x, int y) {
        x = root(x);
        y = root(y);
        parent[x] = y;
    }
};

int main(int argc, char *argv[]) {
    std::string ref_file, runs_file, dump_file;
    int threshold = 200;

    CLI::App args;
    args.add_option("-r", ref_file)->required();
    args.add_option("-l", runs_file)->required();
    args.add_option("-d", dump_file)->required();
    args.add_option("-t", threshold);
    CLI11_PARSE(args, argc, argv);

    core::Dict refs, runs;
    refs.load_file(ref_file);
    printf("loaded \"%s\".\n", ref_file.data());
    runs.load_file(runs_file);
    printf("loaded \"%s\".\n", runs_file.data());

    std::fstream fp(dump_file);

    std::unordered_map<Key, std::vector<core::Vec2d>> bucket;
    std::unordered_map<std::string, std::vector<core::Vec2d>> invs;
    while (fp) {
        std::string data;
        std::getline(fp, data);

        if (data.empty())
            continue;

        std::stringstream stream(data);
        get<std::string>(stream);
        auto name = get<std::string>(stream);

        double x1 = get<int>(stream);
        double y1 = get<int>(stream);

        if (x1 > 0)
            bucket[{name, EndpointType::LEFT}].push_back({x1, y1});

        double x2 = get<int>(stream);
        double y2 = get<int>(stream);
        if (x2 > 0)
            bucket[{name, EndpointType::RIGHT}].push_back({x2, y2});

        auto rate = get<double>(stream);
        if (rate > 0.5)
            invs[name].push_back({x1, x2});
    }

    UnionFind set;
    for (auto &[name, vs] : invs) {
        int n = vs.size();
        set.reset(n);

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if ((vs[i] - vs[j]).len2() < 1e3)
                    set.link(i, j);
            }
        }

        std::vector<core::Vec2d> sum;
        std::vector<int> count;
        sum.resize(n);
        count.resize(n);

        for (int i = 0; i < n; i++) {
            int r = set.root(i);
            sum[r] = sum[r] + vs[i];
            count[r]++;
        }

        for (int i = 0; i < n; i++) {
            if (set.root(i) == i) {
                auto avg = sum[i] / count[i];
                fprintf(stderr,
                    "INV %s %d %d\n", name.data(),
                    int(std::round(avg.x)), int(std::round(avg.y))
                );
            }
        }
    }

    std::unordered_map<Key, std::vector<core::Vec2d>> endpoint;
    for (auto &[key, vs] : bucket) {
        int n = vs.size();
        // fprintf(
        //     stderr, "name=%s, type=%s\n",
        //     key.name.data(), key.type == EndpointType::LEFT ? "LEFT" : "RIGHT"
        // );

        // printf("ListPlot[Style[{{100000,500000},{102000,502000},{104000,504000},{106000,506000},{108000,508000}");
        // for (auto [x, y] : vs) {
        //     printf(",{%.1lf,%.1lf}", x, y);
        // }
        // printf("},RGBColor[0,0,0,0.2],PointSize[0.005]],ImageSize->Full]\n");

        set.reset(n);
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (std::abs(vs[i].x - vs[j].x) <= threshold)
                // if ((vs[i] - vs[j]).len2() <= 9e4)
                    set.link(i, j);
            }
        }

        std::vector<std::vector<core::Vec2d>> list;
        std::vector<int> count;
        list.resize(n);
        count.resize(n);

        for (int i = 0; i < n; i++) {
            int r = set.root(i);
            list[r].push_back(vs[i]);
            count[r]++;
        }

        std::vector<core::Vec2d> ps;
        for (int i = 0; i < n; i++) {
            if (set.root(i) == i) {
                core::Vec2d avg0;
                for (auto p : list[i])
                    avg0 = avg0 + p;
                avg0 = avg0 / count[i];

                auto avg = avg0;
                if (list[i].size() > 1) {
                    auto max_dist = std::numeric_limits<double>::min();
                    int pos = 0;
                    for (int j = 0; j < list[i].size(); j++) {
                        auto dist = std::abs(list[i][j].x - avg0.x);
                        if (dist > max_dist) {
                            max_dist = dist;
                            pos = j;
                        }
                    }

                    avg = core::Vec2d();
                    for (int j = 0; j < list[i].size(); j++) {
                        if (j != pos)
                            avg = avg + list[i][j];
                    }

                    avg = avg / (count[i] - 1);
                }

                ps.push_back(avg);
            }
        }

        std::sort(ps.begin(), ps.end());
        for (auto [fx, fy] : ps) {
            int x = std::round(fx);
            int y = std::round(fy);

            // printf(
            //     "%s %s %d %d\n",
            //     key.name.data(),
            //     key.type == EndpointType::LEFT ? "L" : "R",
            //     x, y
            // );

            endpoint[key].push_back({fx, fy});
        }
    }

    std::unordered_map<std::string, std::vector<core::Range>> ranges;
    for (auto &ref : refs) {
        for (auto lp : endpoint[{ref.name, EndpointType::LEFT}]) {
            for (auto rp : endpoint[{ref.name, EndpointType::RIGHT}]) {
                if (std::abs(lp.x - rp.x) < 50) {
                    int l = (lp.x + rp.x) / 2;
                    int r = l + std::abs(rp.y - lp.y) * 1.2;
                    if (l + 1010 > r)
                        fprintf(stderr, "INS %s %d %d\n", ref.name.data(), l, r);
                } else if (std::abs(lp.x - rp.x) < 1010) {
                    int l = lp.x, r = rp.x;

                    if (std::abs(lp.y - rp.y) < 500) {
                        if (l < r)
                            fprintf(stderr, "DEL %s %d %d\n", ref.name.data(), l, r);
                        else
                            fprintf(stderr, "DUP %s %d %d\n", ref.name.data(), r, l);
                    }

                    if (l < r)
                        ranges[ref.name].push_back({l, r});
                }
            }
        }
    }

    for (int i = 0; i < refs.size(); i++) {
        auto &s1 = refs[i];
        for (auto r1 : ranges[s1.name]) {
            for (int j = i + 1; j < refs.size(); j++) {
                auto &s2 = refs[j];
                for (auto r2 : ranges[s2.name]) {
                    if (std::abs(r1.length() - r2.length()) < 20)
                        fprintf(stderr,
                            "TRA %s %d %d %s %d %d\n",
                            s1.name.data(), r1.begin, r1.end,
                            s2.name.data(), r2.begin, r2.end
                        );
                }
            }
        }
    }

    return 0;
}
