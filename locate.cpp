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

        for (int j = 0; j < runs.size(); j++) {
            if (!core::startswith(runs[j].name, idx))
                continue;
        }
    }

    return 0;
}
