#include "core.hpp"

#include "CLI11.hpp"


struct MetaInfo {
    std::string name;
    std::string target;
    int left;
    int right;
    int loss;
    int reversed;
};

auto load_locate_file(const std::string &path) {
    std::fstream fp(path);

    std::unordered_map<std::string, MetaInfo> meta;
    while (fp) {
        MetaInfo line;
        fp >> line.name >> line.target >> line.left >> line.right >> line.loss >> line.reversed;

        if (!line.name.empty() && !line.target.empty())
            meta[line.name] = line;
    }

    return meta;
}

int main(int argc, char *argv[]) {
    std::string ref_file, locate_file, runs_file, target;

    CLI::App args;
    args.add_option("-r", ref_file)->required();
    args.add_option("-p", locate_file)->required();
    args.add_option("-l", runs_file)->required();
    args.add_option("-t", target)->required();
    CLI11_PARSE(args, argc, argv);

    core::Dict refs, runs;
    refs.load_file(ref_file);
    printf("loaded \"%s\".\n", ref_file.data());
    runs.load_file(runs_file);
    printf("loaded \"%s\".\n", runs_file.data());

    auto meta = load_locate_file(locate_file);
    printf("loaded \"%s\".\n", locate_file.data());

    for (auto &run : runs) {
        auto &info = meta[run.name];

        auto rate = 1.0 - double(info.loss) / run.sequence.size();
        // if (rate >= 0.80)
        //     continue;
        // if (run.name != "S1_54")
        //     continue;
        if (run.name != target)
            continue;

        if (info.reversed)
            run.sequence = core::watson_crick_complement(run.sequence);

        auto &ref = *refs.find(info.target);
        auto s = core::BioSeq(ref.sequence, info.left, info.right + 1);
        auto t = core::BioSeq(run.sequence);

        auto prefix = core::prefix_span(s, t);
        auto suffix = core::suffix_span(s, t);

        printf(
            "%s @%s[%d, %d]: %%=%.3lf, n=%d, m=%d, [%d, %d)-[%d, %d)=%d, [%d, %d)-[%d, %d)=%d\n",
            run.name.data(),
            ref.name.data(),
            info.left, info.right,
            rate, s.size(), t.size(),
            prefix.range1.begin, prefix.range1.end,
            prefix.range2.begin, prefix.range2.end,
            prefix.loss,
            suffix.range1.begin, suffix.range1.end,
            suffix.range2.begin, suffix.range2.end,
            suffix.loss
        );
        printf(
            "(%d, %d, ~%d) (%d, %d, ~%d)\n",
            info.left + prefix.range1.end - 1, prefix.range2.end, prefix.range2.length(),
            info.left + suffix.range1.begin - 1, suffix.range2.begin, suffix.range2.length()
        );
    }

    return 0;
}
