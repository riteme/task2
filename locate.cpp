#include <sstream>

#include "CLI11.hpp"

#include "core.hpp"


auto get_id(int i) -> std::string {
    std::stringstream buffer;
    buffer << 'S' << i;
    return buffer.str();
}

int main(int argc, char *argv[]) {
    std::string ref_path, runs_path;

    CLI::App args;
    args.add_option("-r,--ref", ref_path)->required();
    args.add_option("-l,--long", runs_path)->required();
    CLI11_PARSE(args, argc, argv);

    core::Dict ref, runs;
    ref.load_file(ref_path);
    ref.sort_by_name();
    printf("loaded: \"%s\".\n", ref_path.data());
    runs.load_file(runs_path);
    printf("loaded: \"%s\".\n", runs_path.data());

    for (int i = 0; i < ref.size(); i++) {
        auto idx = get_id(i + 1);
        printf("locating shotguns %s_*...\n", idx.data());

        core::Index index;
        index.append(ref[i].sequence);
        index.build();
        printf("index built for %s.\n", ref[i].name.data());

        double max_match_rate = std::numeric_limits<double>::min();
        double min_match_rate = std::numeric_limits<double>::max();
        for (int j = 0; j < runs.size(); j++) {
            if (!core::startswith(runs[j].name, idx))
                continue;

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
            max_match_rate = std::max(max_match_rate, match_rate);
            min_match_rate = std::min(min_match_rate, match_rate);

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
        }

        fprintf(stderr,
            "%s: min=%.3lf, max=%.3lf\n",
            ref[i].name.data(),
            min_match_rate * 100,
            max_match_rate * 100
        );
    }

    return 0;
}
