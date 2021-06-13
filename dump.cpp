#include "CLI11.hpp"

#include "core.hpp"
#include "rash/pool.hpp"


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
    double max_rate = 0.84;
    int n_workers = 1;

    CLI::App args;
    args.add_option("-r", ref_file)->required();
    args.add_option("-p", locate_file)->required();
    args.add_option("-l", runs_file)->required();
    args.add_option("-t", target);
    args.add_option("-m", max_rate);
    args.add_option("-j", n_workers);
    CLI11_PARSE(args, argc, argv);

    core::Dict refs, runs;
    refs.load_file(ref_file);
    printf("loaded \"%s\".\n", ref_file.data());
    runs.load_file(runs_file);
    printf("loaded \"%s\".\n", runs_file.data());

    auto meta = load_locate_file(locate_file);
    printf("loaded \"%s\".\n", locate_file.data());

    ThreadPool pool(n_workers);
    std::vector<std::future<void>> futures;
    futures.reserve(runs.size());
    for (auto &run : runs) {
        auto future = pool.run([&]() {
            auto &info = meta[run.name];

            auto rate = 1.0 - double(info.loss) / run.sequence.size();
            if (rate > max_rate ||
                (!target.empty() && run.name != target))
                return;

            if (info.reversed)
                run.sequence = core::watson_crick_complement(run.sequence);

            auto &ref = *refs.find(info.target);
            auto s = core::BioSeq(ref.sequence, info.left, info.right + 1);
            auto t = core::BioSeq(run.sequence);

            auto prefix = core::prefix_span(s, t);
            auto suffix = core::suffix_span(s, t);

            bool contained = true;
            core::Vec2i front, back;

            if (100 < prefix.range2.length() && prefix.range2.length() < run.sequence.size() - 100) {
                front = core::Vec2i(
                    info.left + prefix.range1.end - 1,
                    info.left + prefix.range2.end - 1
                );
            } else
                contained = false;

            if (100 < suffix.range2.length() && suffix.range2.length() < run.sequence.size() - 100) {
                back = core::Vec2i(
                    info.left + suffix.range1.begin - 2,
                    info.left + suffix.range2.begin - 2
                );
            } else
                contained = false;

            auto inv_match_rate = -1.0;
            int l1 = front.x, r1 = back.x;
            int l2 = prefix.range2.end, r2 = suffix.range2.begin - 1;
            if (contained && 0 < l1 && l1 <= r1 && 0 < l2 && l2 <= r2) {
                auto reversed = core::watson_crick_complement(
                    run.sequence.substr(l2 - 1, r2 - l2 + 1)
                );

                int loss = core::full_align(
                    core::BioSeq(ref.sequence, front.x, back.x + 1),
                    reversed
                );

                inv_match_rate = 1 - 2.0 * loss / (r1 - l1 + 1 + r2 - l2 + 1);
            }

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
            fprintf(stderr,
                "%s %s %d %d %d %d %.16lf\n",
                run.name.data(),
                ref.name.data(),
                front.x, front.y,
                back.x, back.y,
                inv_match_rate
            );
        });

        futures.push_back(std::move(future));
    }

    for (auto &f : futures) {
        f.get();
    }

    return 0;
}
