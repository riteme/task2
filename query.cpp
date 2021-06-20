#include <sstream>
#include <fstream>

#include "CLI11.hpp"


int main(int argc, char *argv[]) {
    std::string locate_file;
    int left, right;

    CLI::App args;
    args.add_option("-p", locate_file)->required();
    args.add_option("-l", left)->required();
    args.add_option("-r", right)->required();
    CLI11_PARSE(args, argc, argv);

    if (left > right) {
        fprintf(stderr, "left=%d should not be greater than right=%d.\n", left, right);
        return -1;
    }

    std::fstream fp(locate_file);
    while (fp) {
        std::string line;
        std::getline(fp, line);

        std::stringstream s(line);

        std::string name, target;
        int l, r;
        s >> name >> target >> l >> r;

        bool intersected =
            (left <= l && l <= right) ||
            (left <= r && r <= right) ||
            (l <= left && left <= r) ||
            (l <= right && right <= r);

        if (intersected) {
            struct Item {
                char id;
                int position;

                bool operator<(const Item &rhs) const {
                    return position < rhs.position;
                }
            };

            Item a[] = {
                {'A', left}, {'A', right},
                {'B', l}, {'B', r}
            };
            std::sort(std::begin(a), std::end(a));

            printf("%10s @%-16s [%d, %d]: ", name.data(), target.data(), l, r);
            for (int i = 0; i < 4; i++) {
                if (i > 0)
                    printf(" --%d-- ", a[i].position - a[i - 1].position);
                printf("%c", a[i].id);
            }
            puts("");
        }
    }

    return 0;
}
