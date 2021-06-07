#include "common.hpp"


namespace core {

bool startswith(const std::string &target, const std::string &pattern) {
    if (target.size() < pattern.size())
        return false;

    for (size_t i = 0; i < pattern.size(); i++) {
        if (target[i] != pattern[i])
            return false;
    }

    return true;
}

auto watson_crick_complement(const std::string &s) -> std::string {
    std::string t(s.rbegin(), s.rend());

    for (size_t i = 0; i < t.size(); i++) {
        switch (t[i]) {
            case 'A': t[i] = 'T'; break;
            case 'T': t[i] = 'A'; break;
            case 'C': t[i] = 'G'; break;
            case 'G': t[i] = 'C'; break;
            default: t[i] = 'N';
        }
    }

    return t;
}

}
