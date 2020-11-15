cmake_minimum_required(VERSION 3.5)
project(aog_lib)

set(CMAKE_CXX_STANDARD 11)

find_package(Boost COMPONENTS filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIR})

# add_subdirectory(Testcases)
# add_subdirectory(Faces)
add_subdirectory(AOF_Test)
