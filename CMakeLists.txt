cmake_minimum_required(VERSION 3.20)

project(task2)

# set(CMAKE_BUILD_TYPE Debug)
# set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=undefined,address")
# set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=undefined,address")

set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION True)

include(CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

CHECK_CXX_COMPILER_FLAG("-mtune=native" COMPILER_SUPPORTS_MTUNE_NATIVE)
if(COMPILER_SUPPORTS_MTUNE_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

find_library(OpenGL GL REQUIRED)
find_library(SDL2 SDL2 REQUIRED)
find_library(GLEW GLEW REQUIRED)
find_library(pthread pthread REQUIRED)

include_directories(include)
include_directories(thirdparty)

file(GLOB core_sources CONFIGURE_DEPENDS source/*.cpp)
add_library(core STATIC ${core_sources})

add_subdirectory(thirdparty/imgui)
add_subdirectory(thirdparty/nanovg)
add_subdirectory(thirdparty/rash)

add_executable(align align.cpp)
add_executable(locate locate.cpp)
add_executable(locate-demo locate-demo.cpp)
add_executable(dump dump.cpp)
add_executable(aggregate aggregate.cpp)
add_executable(analyze analyze.cpp)
add_executable(query query.cpp)

set(cxx_options
    -Wall
    -Wextra
    -Wno-char-subscripts
    -Wno-sign-compare
)

target_compile_options(core PRIVATE ${cxx_options})
target_compile_options(align PRIVATE ${cxx_options})
target_compile_options(locate PRIVATE ${cxx_options})
target_compile_options(locate-demo PRIVATE ${cxx_options})
target_compile_options(dump PRIVATE ${cxx_options})
target_compile_options(aggregate PRIVATE ${cxx_options})
target_compile_options(analyze PRIVATE ${cxx_options})
target_compile_options(query PRIVATE ${cxx_options})

target_link_libraries(align core)
target_link_libraries(locate core rash pthread)
target_link_libraries(locate-demo core rash imgui nanovg OpenGL GLEW SDL2 pthread)
target_link_libraries(dump core rash pthread)
target_link_libraries(aggregate core)
target_link_libraries(analyze core)
