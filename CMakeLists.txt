set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
cmake_minimum_required(VERSION 3.25)
project(tree LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(tree
    src/main.cpp
    src/galois_field.cpp
    src/algorithm.cpp
    tests/test.cpp
)
