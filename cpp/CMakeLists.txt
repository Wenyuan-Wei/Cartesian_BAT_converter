cmake_minimum_required(VERSION 3.15)
project(cbc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(cbc STATIC
    core/vec3_ops.cpp
    core/geometry.cpp
    core/structure.cpp
    src/pdb_reader.cpp
)

target_include_directories(cbc
  PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/core
)

add_executable(test test.cpp)
target_link_libraries(test PRIVATE cbc)