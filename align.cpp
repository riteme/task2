#include <string>
#include <fstream>

#include "thirdparty/CLI11.hpp"

#include "core.hpp"


void print_set(const std::vector<int> &set) {
    putchar('[');
    for (size_t i = 0; i < set.size(); i++) {
        printf("%d", set[i]);
        if (i + 1 < set.size())
            putchar(' ');
    }
    puts("]");
}

int main(int argc, char *argv[]) {
    CLI::App args;

    std::string ref_fasta, long_fasta;
    args.add_option("-r,--ref", ref_fasta)->required();
    args.add_option("-l,--long", long_fasta)->required();

    CLI11_PARSE(args, argc, argv);

    core::Dict ref, runs;
    ref.load_file(ref_fasta);
    printf("read %zu string(s) from \"%s\".\n", ref.size(), ref_fasta.data());
    runs.load_file(long_fasta);
    printf("read %zu string(s) from \"%s\".\n", runs.size(), long_fasta.data());

    core::Index index;
    for (auto &e : ref) {
        index.append(e.sequence);
        index.append(core::CMAP['N']);
    }
    index.build();
    printf("index built.\n");

    while (true) {
        int idx, left, right;
        printf("select run id [0-%zu]: ", runs.size() - 1);
        if (scanf("%d", &idx) == EOF)
            break;

        printf("select range [l, r] [1-%zu]: ", runs[idx].sequence.size());
        if (scanf("%d %d", &left, &right) == EOF)
            break;

        std::string pattern = runs[idx].sequence;
        if (right < pattern.size())
            pattern.erase(right + 1);
        if (left > 1)
            pattern.erase(0, left - 1);

        // std::string pattern = "TCGAATCGCCAAGGTAAGGTGTAAAATGGTGATAGCGTTCCCCATCAATCTCGTAGCCAGTGTAGTCCATAAGCCTCTCTAAGTAACCCCCTTGATAAGCTTCA";

        // auto t = index.locate(pattern);
        // auto rp = index.rpset(t);
        // printf("t.id=%d, t.len=%d\n", t.id, t.len);
        // print_set(rp);

        auto r = index.align(pattern);

        if (r.token.id > 1) {
            auto rp = index.rpset(r.token);
            print_set(rp);
        } else {
            puts("cancelled.");
        }
    }

    return 0;
}
