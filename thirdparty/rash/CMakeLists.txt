cmake_minimum_required(VERSION 3.20)

file(GLOB rash_sources CONFIGURE_DEPENDS "*.cpp")

add_library(rash STATIC ${rash_sources})
