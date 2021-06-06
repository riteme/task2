#include "core.hpp"
#include "thirdparty/rash/rash.hpp"


constexpr auto NVG_FONT_SIZE = 25.0f;

int main() {
    Application app("Alignment Visualization");

    bool show_demo = false;
    char ref_path[1024] = "";
    char runs_path[1024] = "";
    core::Dict ref, runs;
    int run_id = 0;

    app.run([&] {
        app.gl_clear(1.0f, 1.0f, 1.0f);

        auto vg = app.get_nanovg_context();
        nvgBeginFrame(vg, WINDOW_WIDTH, WINDOW_HEIGHT, 1.0f);

        for (int i = 0; i < ref.size(); i++) {
            float top = i * 100.0f;

            nvgBeginPath(vg);
            app.nvg_text_mode(nvgRGB(0, 0, 0), NVG_FONT_SIZE, NVG_ALIGN_LEFT | NVG_ALIGN_BOTTOM);
            nvgText(vg, 0.0f, top + 30.0f, ref[i].name.data(), NULL);
            nvgFill(vg);

            nvgBeginPath(vg);
            nvgRect(vg, 0, top + 31.0f, WINDOW_WIDTH, 20.0f);
            nvgFillColor(vg, nvgRGB(33, 150, 243));
            nvgFill(vg);
        }

        nvgEndFrame(vg);

        ImGui::Begin("Workflow");

        if (ImGui::Button("Show Demo"))
            show_demo = !show_demo;
        if (show_demo)
            ImGui::ShowDemoWindow(&show_demo);

        if (ImGui::Button("Open...##1"))
            zenity_file_selection_fgets("Select reference file", ref_path, 1024);
        ImGui::SameLine();
        ImGui::InputText("ref", ref_path, 1024);

        if (ImGui::Button("Open...##2"))
            zenity_file_selection_fgets("Select run file", runs_path, 1024);
        ImGui::SameLine();
        ImGui::InputText("runs", runs_path, 1024);

        if (ImGui::Button("Load")) {
            puts(ref_path);
            ref.load_file(ref_path);
            puts(runs_path);
            runs.load_file(runs_path);
        }

        int int_zero = 0;
        int run_id_max = runs.size() - 1;
        bool loaded = ref.size() > 0 && runs.size() > 0;
        if (loaded) {
            ImGui::DragScalar("run id", ImGuiDataType_S32, &run_id, 1.0f, &int_zero, &run_id_max);
            ImGui::SameLine();
            ImGui::Text("[0-%d]", run_id_max);

            ImGui::Button("Run");

            ImGui::Text("run: %s", runs[run_id].name.data());
        }

        ImGui::End();
    });

    return 0;
}
