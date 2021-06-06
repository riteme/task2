#include <map>
#include <mutex>
#include <atomic>
#include <vector>
#include <thread>
#include <fstream>

#include "core.hpp"
#include "thirdparty/rash/rash.hpp"


constexpr int K = 20;
constexpr auto NVG_FONT_SIZE = 25.0f;

int main() {
    Application app("Alignment Visualization");

    bool show_demo = false;
    char ref_path[1024] = "";
    char runs_path[1024] = "";
    char bed_path[1024] = "";
    core::Dict ref, runs;
    int run_id = 0;
    std::mutex lock;
    std::vector<std::vector<int>> sv;
    std::vector<std::map<Vec2i, int>> locations;
    std::vector<int> positions;
    std::vector<core::Index *> indexes;
    float marker_alpha = 3.0f;

    auto reset_markers = [&] {
        lock.lock();
        locations.clear();
        locations.resize(ref.size());
        lock.unlock();

        positions.clear();
        positions.resize(ref.size() * 2);
    };

    auto load_files = [&] {
        puts(ref_path);
        ref.load_file(ref_path);
        puts(runs_path);
        runs.load_file(runs_path);

        sv.clear();
        sv.resize(ref.size());
        std::fstream fp(bed_path);

        auto read_one = [&] {
            std::string name, left_literal, right_literal;
            fp >> name >> left_literal >> right_literal;
            printf("%s %s %s\n", name.data(), left_literal.data(), right_literal.data());

            int left = std::stoi(left_literal);
            int right = std::stoi(right_literal);
            for (int i = 0; i < ref.size(); i++) {
                if (ref[i].name == name) {
                    int middle = (left + right) / 2;
                    float hv = float(middle) / ref[i].sequence.size();
                    sv[i].push_back(hv * WINDOW_WIDTH);
                    break;
                }
            }
        };

        while (fp) {
            std::string op;
            fp >> op;
            if (op.empty())
                break;

            read_one();
            if (op == "TRA")
                read_one();
        }

        reset_markers();

        for (auto i : indexes) {
            delete i;
        }
        indexes.clear();

        for (auto &e : ref) {
            auto i = new core::Index;
            i->append(e.sequence);
            i->build();
            indexes.push_back(i);
        }
    };

    std::atomic<bool> locate_disabled;

    auto locate = [&](int id, int ref_id, int y, bool rev) {
        std::string seq;
        seq.resize(K);

        auto index = indexes[ref_id];
        auto run = runs[run_id];

        if (rev) {
            auto &s = run.sequence;
            std::reverse(s.begin(), s.end());
            for (int i = 0; i < s.size(); i++) {
                if (s[i] == 'A')
                    s[i] = 'T';
                else if (s[i] == 'T')
                    s[i] = 'A';
                else if (s[i] == 'C')
                    s[i] = 'G';
                else if (s[i] == 'G')
                    s[i] = 'C';
            }
        }

        for (int i = 0; i + K <= run.sequence.size(); i++) {
            if (locate_disabled)
                return;

            for (int j = 0; j < K; j++) {
                seq[j] = run.sequence[i + j];
            }

            auto alignment = index->align(seq);
            auto set = index->rpset(alignment.token);

            lock.lock();
            float len = ref[ref_id].sequence.size();
            for (int rp : set) {
                float hv = float(std::max(0, rp - K / 2)) / len;
                int cp = hv * WINDOW_WIDTH;
                locations[ref_id][{cp, y}] += 1;
            }
            lock.unlock();

            positions[id] = i + 1;
        }
    };

    app.run([&] {
        int mouse_x = app.get_mouse_position().x;
        float mouse_hv = float(mouse_x) / WINDOW_WIDTH;

        app.gl_clear(1.0f, 1.0f, 1.0f);

        auto vg = app.get_nanovg_context();
        nvgBeginFrame(vg, WINDOW_WIDTH, WINDOW_HEIGHT, 1.0f);

        for (int i = 0; i < ref.size(); i++) {
            float top = i * 150.0f;

            nvgBeginPath(vg);
            app.nvg_text_mode(nvgRGB(0, 0, 0), NVG_FONT_SIZE, NVG_ALIGN_LEFT | NVG_ALIGN_BOTTOM);
            nvgText(vg, 0.0f, top + 30.0f, ref[i].name.data(), NULL);
            nvgFill(vg);

            nvgBeginPath(vg);
            nvgRect(vg, 0, top + 31.0f, WINDOW_WIDTH, 20.0f);
            nvgFillColor(vg, nvgRGB(33, 150, 243));
            nvgFill(vg);

            for (int x : sv[i]) {
                nvgBeginPath(vg);
                nvgRect(vg, x - 5.0f, top + 33.0f, 10.0f, 16.0f);
                nvgFillColor(vg, nvgRGBA(255, 0, 0, 100));
                nvgFill(vg);
            }

            lock.lock();
            for (auto &e : locations[i]) {
                auto [x, y] = e.first;
                int count = e.second;

                nvgBeginPath(vg);
                nvgRect(vg, x - 5.0f, top + 55.0f + y, 10.0f, 10.0f);
                nvgFillColor(vg, nvgRGBA(0, 0, 0, std::min(255, int(marker_alpha * count))));
                nvgFill(vg);
            }
            lock.unlock();
        }

        nvgBeginPath(vg);
        nvgMoveTo(vg, mouse_x, 0.0f);
        nvgLineTo(vg, mouse_x, WINDOW_HEIGHT);
        nvgStrokeColor(vg, nvgRGB(255, 152, 0));
        nvgStroke(vg);

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

        if (ImGui::Button("Open...##3"))
            zenity_file_selection_fgets("Select sv.bed", bed_path, 1024);
        ImGui::SameLine();
        ImGui::InputText("sv.bed", bed_path, 1024);

        if (ImGui::Button("Load"))
            load_files();

        int int_zero = 0;
        int run_id_max = runs.size() - 1;
        bool loaded = ref.size() > 0 && runs.size() > 0;
        if (loaded) {
            ImGui::InputInt("run id", &run_id);
            if (run_id < 0)
                run_id = 0;
            if (run_id > run_id_max)
                run_id = run_id_max;
            ImGui::SameLine();
            ImGui::Text("[0-%d]", run_id_max);

            float float0 = 0.0f, float255 = 255.0f;
            ImGui::DragScalar("marker alpha", ImGuiDataType_Float, &marker_alpha, 0.02f, &float0, &float255);

            if (ImGui::Button("Run")) {
                reset_markers();

                for (int i = 0; i < ref.size(); i++) {
                    std::thread t1(locate, 2 * i + 0, i, 0, false);
                    std::thread t2(locate, 2 * i + 1, i, 20, true);
                    t1.detach();
                    t2.detach();
                }
            }
            ImGui::SameLine();
            if (ImGui::Button("Reset"))
                reset_markers();
            ImGui::SameLine();
            if (ImGui::Button("Enable"))
                locate_disabled = false;
            ImGui::SameLine();
            if (ImGui::Button("Disable"))
                locate_disabled = true;
            ImGui::SameLine();
            ImGui::Text("%s", locate_disabled ? "disabled" : "enabled");

            ImGui::Text("run: %s, len=%zu, K=%d",
                runs[run_id].name.data(), runs[run_id].sequence.size(), K
            );

            for (int i = 0; i < positions.size(); i++) {
                ImGui::Text("[%d]: %d", i, positions[i]);
            }
        }

        ImGui::Text("mouse.x = %d (%.3f)", mouse_x, mouse_hv);
        for (int i = 0; i < ref.size(); i++) {
            ImGui::Text("%s: %.1f", ref[i].name.data(), ref[i].sequence.size() * mouse_hv);
        }

        ImGui::End();
    });

    return 0;
}
