#include "rash.hpp"

#include "imgui/imgui_impl_sdl.h"
#include "imgui/imgui_impl_opengl3.h"

#define NANOVG_GL3_IMPLEMENTATION
#include "nanovg/nanovg_gl.h"


Application::Application(const std::string &title) {
    // SDL: window and OpenGL context.
    if (SDL_Init(SDL_INIT_VIDEO) != 0)
        abort_with_sdl_error();

    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

    _window = SDL_CreateWindow(
        title.data(),
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        WINDOW_WIDTH, WINDOW_HEIGHT,
        SDL_WINDOW_OPENGL | SDL_WINDOW_ALLOW_HIGHDPI
    );
    if (_window == NULL)
        abort_with_sdl_error();

    _context = SDL_GL_CreateContext(_window);
    if (_context == NULL)
        abort_with_sdl_error();

    SDL_GL_MakeCurrent(_window, _context);

    // GLEW
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK)
        abort_with_message("failed to init GLEW");

    _vg = nvgCreateGL3(NVG_ANTIALIAS | NVG_STENCIL_STROKES);
    if (_vg == NULL)
        abort_with_message("failed to create nanovg context");

    ImGui::CreateContext();
    ImGui::StyleColorsLight();

    load_fonts();

    ImGui_ImplSDL2_InitForOpenGL(_window, _context);
    ImGui_ImplOpenGL3_Init();

    SDL_GL_SetSwapInterval(1);
}

Application::~Application() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();
    nvgDeleteGL3(_vg);
    SDL_GL_DeleteContext(_context);
    SDL_DestroyWindow(_window);
    SDL_Quit();
}

void Application::load_fonts() {
    ImGuiIO &io = ImGui::GetIO();

    bool first = true;
    ImFontConfig font_config;

    for (auto font : FONT_FILES) {
        font_config.MergeMode = !first;
        if (first) {
            _default_font = nvgCreateFont(_vg, "default", font.path);
            if (_default_font == -1)
                abort_with_message("failed to load default font: %s", font.path);
            first = false;
        } else {
            int new_font = nvgCreateFont(_vg, font.path, font.path);
            if (new_font == -1)
                abort_with_message("failed to load fallback font: %s", font.path);
            if (!nvgAddFallbackFontId(_vg, _default_font, new_font))
                abort_with_message("failed to add fallback font to default font");
        }

        if (font.is_cjk)
            io.Fonts->AddFontFromFileTTF(
                font.path, font.size, &font_config,
                io.Fonts->GetGlyphRangesChineseSimplifiedCommon()
            );
        else
            io.Fonts->AddFontFromFileTTF(
                font.path, font.size, &font_config
            );
    }

    io.Fonts->Build();
}

void Application::run(const std::function<void()> &main) {
    while (!_stopped) {
        poll_events();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame(_window);
        ImGui::NewFrame();

        setup_gl_viewport();
        main();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        SDL_GL_SwapWindow(_window);
    }
}

void Application::poll_events() {
    SDL_Event e;
    while (SDL_PollEvent(&e)) {
        ImGui_ImplSDL2_ProcessEvent(&e);
        if (e.type == SDL_QUIT ||
            (e.type == SDL_WINDOWEVENT &&
                e.window.event == SDL_WINDOWEVENT_CLOSE &&
                e.window.windowID == SDL_GetWindowID(_window)))
            _stopped = true;
    }
}

void Application::setup_gl_viewport() {
    ImGuiIO &io = ImGui::GetIO();
    glViewport(0, 0, io.DisplaySize.x, io.DisplaySize.y);
    glMatrixMode(GL_PROJECTION_MATRIX);
    glLoadIdentity();
    glOrtho(0, io.DisplaySize.x, io.DisplaySize.y, 0, 0, 1);
}

void Application::gl_clear(float r, float g, float b) {
    glClearColor(r, g, b, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
}

auto Application::get_nanovg_context() -> NVGcontext * {
    return _vg;
}

void Application::nvg_text_mode(NVGcolor color, float size, int align) {
    nvgFillColor(_vg, color);
    nvgFontSize(_vg, size);
    nvgFontFaceId(_vg, _default_font);
    nvgTextAlign(_vg, align);
}

auto Application::get_mouse_position() -> Vec2i {
    int x, y;
    SDL_GetMouseState(&x, &y);
    return {x, y};
}
