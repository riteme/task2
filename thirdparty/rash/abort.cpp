#include "rash.hpp"

#include <cstdlib>


void abort_with_message(const char *message, ...) {
    va_list args;
    va_start(args, message);
    vfprintf(stderr, message, args);
    va_end(args);
    abort();
}

void abort_with_errno() {
    fprintf(stderr, "ERRNO: %s\n", strerror(errno));
    abort();
}

void abort_with_sdl_error() {
    fprintf(stderr, "SDL error: %s\n", SDL_GetError());
    abort();
}
