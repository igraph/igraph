cmake_minimum_required(VERSION 3.21)
project(infomap)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Wshadow")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP)

find_program(CCACHE ccache)
if(CCACHE)
    message("-- Found Ccache: TRUE")
    set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
else()
    message("-- Found Ccache: FALSE")
endif()

file(GLOB_RECURSE infomap_src CONFIGURE_DEPENDS "src/*.h" "src/*.cpp")

add_executable(infomap ${infomap_src})

if(OpenMP_CXX_FOUND)
    target_link_libraries(infomap PUBLIC OpenMP::OpenMP_CXX)
endif()
