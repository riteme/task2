#pragma once

#include <string>
#include <functional>

#include "imgui/imgui.h"

#include <SDL2/SDL.h>

#include <GL/glew.h>
#include <SDL2/SDL_opengl.h>

#include "nanovg/nanovg.h"


constexpr auto WINDOW_WIDTH = 1920;
constexpr auto WINDOW_HEIGHT = 1080;

struct FontInfo {
    const char *path;
    bool is_cjk;
    float size;
};

constexpr FontInfo FONT_FILES[] = {
    {"../asset/FiraMono-Regular.otf", false, 30.0f},
    {"../asset/SourceHanSansCN-Regular.otf", true, 37.0f},
};

void abort_with_message(const char *message, ...);
void abort_with_errno();
void abort_with_sdl_error();

auto zenity_file_selection(const std::string &title) -> std::string;
void zenity_file_selection_fgets(const std::string &title, char *buffer, int n);

struct Vec2i {
    int x, y;

    bool operator<(const Vec2i &rhs) const {
        return std::tie(x, y) < std::tie(rhs.x, rhs.y);
    }
};

class Application {
public:
    Application(const std::string &title);
    ~Application();

    Application(const Application &) = delete;
    Application(Application &&) = delete;
    auto operator=(const Application &) = delete;
    auto operator=(Application &&) = delete;

    void load_fonts();
    void run(const std::function<void()> &main);
    void poll_events();
    void setup_gl_viewport();
    void gl_clear(float r, float g, float b);
    auto get_nanovg_context() -> NVGcontext *;
    void nvg_text_mode(NVGcolor color, float size, int align);
    auto get_mouse_position() -> Vec2i;

private:
    bool _stopped = false;
    SDL_Window *_window = NULL;
    SDL_GLContext _context = NULL;
    NVGcontext *_vg = NULL;
    int _default_font = 0;
};

#include "pool.hpp"
