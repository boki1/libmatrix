cmake_minimum_required(VERSION 3.17)
project(libmatrix)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(libmatrix
        main.cpp
        matrix.hpp
        test.cpp test.hpp matrix_impl.hpp)
