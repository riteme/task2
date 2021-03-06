#pragma once

#include <vector>
#include <string>
#include <unordered_map>


namespace core {

struct DictEntry {
    DictEntry() = default;
    DictEntry(const DictEntry &) = default;
    DictEntry(DictEntry &&) = default;
    DictEntry(const std::string &_name, const std::string &_sqeuence)
        : name(_name), sequence(_sqeuence) {}

    DictEntry &operator=(const DictEntry &) = default;
    DictEntry &operator=(DictEntry &&) = default;

    std::string name;
    std::string sequence;
};

// load a single FASTA file into an ordered list of dict entries.
class Dict {
public:
    void load_file(const std::string &path);
    void sort_by_name();
    void build_index();
    auto find(const std::string &name) const -> DictEntry *;

    auto &operator[](size_t i) {
        return _entries[i];
    }

    auto &operator[](size_t i) const {
        return _entries[i];
    }

    auto begin() {
        return _entries.begin();
    }

    auto begin() const {
        return _entries.begin();
    }

    auto end() {
        return _entries.end();
    }

    auto end() const {
        return _entries.end();
    }

    auto size() const {
        return _entries.size();
    }

private:
    auto _find_impl(const std::string &name) -> DictEntry *;

    mutable std::vector<DictEntry> _entries;
    std::unordered_map<std::string, DictEntry *> _index;
};

}
