#include <sstream>

#include "CLI11.hpp"

#include "core.hpp"
#include "rash/rash.hpp"


auto get_id(int i) -> std::string {
    std::stringstream buffer;
    buffer << 'S' << i;
    return buffer.str();
}

int main(int argc, char *argv[]) {
    int n_workers = 1;
    std::string ref_path, runs_path;

    CLI::App args;
    args.add_option("-r,--ref", ref_path)->required();
    args.add_option("-l,--long", runs_path)->required();
    args.add_option("-j,--jobs", n_workers);
    CLI11_PARSE(args, argc, argv);

    core::Dict ref, runs;
    ref.load_file(ref_path);
    ref.sort_by_name();
    printf("loaded: \"%s\".\n", ref_path.data());
    runs.load_file(runs_path);
    printf("loaded: \"%s\".\n", runs_path.data());

    ThreadPool pool(n_workers);

    for (int i = 0; i < ref.size(); i++) {
        auto idx = get_id(i + 1);
        printf("locating shotguns %s_*...\n", idx.data());

        core::Index index;
        index.append(ref[i].sequence);
        index.build();
        printf("index built for %s.\n", ref[i].name.data());

        std::vector<std::future<void>> futures;
        for (int j = 0; j < runs.size(); j++) {
            if (!core::startswith(runs[j].name, idx))
                continue;

            auto future = pool.run([&ref, &runs, &index, i, j] {
                auto &t = runs[j].sequence;
                auto location = index.fuzzy_locate(t);

                auto s = core::BioSeq(ref[i].sequence, location.left, location.right + 1);

                core::Alignment result;
                if (location.reversed) {
                    auto q = core::watson_crick_complement(t);
                    result = core::local_align(s, q);
                } else
                    result = core::local_align(s, t);

                int left = result.range1.begin + location.left - 1;
                int right = result.range1.end - 1 + location.left - 1;
                int length = result.range1.end - result.range1.begin;
                double match_rate = double(t.size() - result.loss) / t.size();

                printf(
                    "%s @%s: [%d, %d], loss=%d (%.3lf%%), ratio=%.3lf, rev=%d\n",
                    runs[j].name.data(),
                    ref[i].name.data(),
                    left, right,
                    result.loss,
                    match_rate * 100,
                    double(length) / runs[j].sequence.size(),
                    location.reversed
                );
                fprintf(stderr,
                    "%s %s %d %d %d %d\n",
                    runs[j].name.data(),
                    ref[i].name.data(),
                    left, right,
                    result.loss,
                    location.reversed
                );
            });

            futures.push_back(std::move(future));
        }

        for (auto &f : futures) {
            f.get();
        }

        printf("%s_*: %s completed.\n", idx.data(), ref[i].name.data());
    }

    return 0;
}
