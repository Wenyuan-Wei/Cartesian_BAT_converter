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
    src/structure/bfs_tree.cpp
    src/structure/protein_root.cpp
    src/io/pdb_reader.cpp
    src/bat/bat_coords.cpp
    src/bat/bat_reconstruct.cpp
    src/io/bat_writer.cpp
)

target_include_directories(cbc
  PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/cartesian
)

add_executable(test test.cpp)
target_link_libraries(test PRIVATE cbc)

add_executable(test_bonds test_bonds.cpp)
target_link_libraries(test_bonds PRIVATE cbc)

add_executable(test_bfs_tree test_bfs_tree.cpp)
target_link_libraries(test_bfs_tree PRIVATE cbc)

add_executable(test_backbone_frame test_backbone_frame.cpp)
target_link_libraries(test_backbone_frame PRIVATE cbc)

add_executable(test_bat_coords test_bat_coords.cpp)
target_link_libraries(test_bat_coords PRIVATE cbc)

add_executable(test_bat_roundtrip test_bat_roundtrip.cpp)
target_link_libraries(test_bat_roundtrip PRIVATE cbc)

add_executable(test_bat_writer test_bat_writer.cpp)
target_link_libraries(test_bat_writer PRIVATE cbc)

add_executable(bat_convert main_bat.cpp)
target_link_libraries(bat_convert PRIVATE cbc)
