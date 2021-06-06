#include "rash.hpp"

#include <sstream>


auto zenity_file_selection(const std::string &title) -> std::string  {
    std::stringstream command;
    command << "zenity --file-selection --title='" << title << "'";
    auto fp = popen(command.str().data(), "r");
    if (fp == NULL)
        abort_with_errno();

    char buffer[512];
    std::stringstream path;
    while (fgets(buffer, sizeof(buffer), fp) != NULL) {
        path << buffer;
    }

    // remove tailing newline.
    auto data = path.str();
    for (size_t i = 0; i < data.size(); i++) {
        if (data[i] == '\n') {
            data.erase(i);
            break;
        }
    }

    return data;
}

void zenity_file_selection_fgets(const std::string &title, char *buffer, int n) {
    std::stringstream command;
    command << "zenity --file-selection --title='" << title << "'";
    auto fp = popen(command.str().data(), "r");
    if (fp == NULL)
        abort_with_errno();

    fgets(buffer, n, fp);

    // remove tailing newline.
    for (int i = 0; i < n; i++) {
        if (buffer[i] == '\n') {
            buffer[i] = 0;
            break;
        }
    }
}
