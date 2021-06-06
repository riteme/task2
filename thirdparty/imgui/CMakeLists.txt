cmake_minimum_required(VERSION 3.20)

file(GLOB imgui_sources CONFIGURE_DEPENDS "*.cpp")

add_library(imgui STATIC ${imgui_sources})