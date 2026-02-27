cmake_minimum_required(VERSION 3.15)
project(cbc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(cbc STATIC
    src/cartesian/vec3_ops.cpp
    src/cartesian/geometry.cpp
    src/structure/structure.cpp
    src/structure/bond_inference.cpp
    src/structure/structure_graph.cpp
    src/io/pdb_reader.cpp
)

target_include_directories(cbc
  PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/cartesian
)

add_executable(test test.cpp)
target_link_libraries(test PRIVATE cbc)

add_executable(test_bonds test_bonds.cpp)
target_link_libraries(test_bonds PRIVATE cbc)