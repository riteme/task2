#include <cctype>
#include <cassert>

#include <fstream>
#include <algorithm>

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
    _index.clear();

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

void Dict::sort_by_name() {
    std::sort(
        _entries.begin(), _entries.end(),
        [](const DictEntry &u, const DictEntry &v) {
            return u.name < v.name;
        }
    );
}

void Dict::build_index() {
    _index.clear();
    for (auto &e : _entries) {
        _index[e.name] = &e;
    }
}

auto Dict::find(const std::string &name) const -> DictEntry * {
    if (_index.empty()) {
        for (auto &e : _entries) {
            if (e.name == name)
                return &e;
        }

        return nullptr;
    } else {
        auto it = _index.find(name);
        if (it == _index.end())
            return nullptr;
        else
            return it->second;
    }
}

}
