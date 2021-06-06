#include <cctype>
#include <cassert>

#include <fstream>

#include "dict.hpp"


namespace {

void trim(std::string &s) {
    size_t i;

    for (i = 0; i < s.size(); i++) {
        if (isalnum(s[i]))
            break;
    }

    if (i > 0)
        s.erase(0, i);

    if (s.size() > 0) {
        assert(isalnum(s[0]));

        for (i = s.size() - 1; ; i--) {
            if (isalnum(s[i]))
                break;
        }

        if (i + 1 < s.size())
            s.erase(i + 1);
    }
}

}

namespace core {

void Dict::load_file(const std::string &path) {
    _entries.clear();

    std::fstream fp(path);
    while (fp) {
        std::string name, sequence;

        do {
            std::getline(fp, name);
        } while (name.empty() && fp);

        trim(name);

        if (name.empty())
            return;

        do {
            assert(fp);
            std::getline(fp, sequence);
        } while (sequence.empty());

        trim(sequence);
        _entries.emplace_back(name, sequence);
    }
}

}
