cmake_minimum_required(VERSION 3.20)

file(GLOB nanovg_sources CONFIGURE_DEPENDS "*.c")

add_library(nanovg STATIC ${nanovg_sources})
