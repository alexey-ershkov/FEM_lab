cmake_minimum_required(VERSION 3.21)
project(mke)

set(CMAKE_CXX_STANDARD 17)

add_executable(mke main.cpp Matrix.cpp Matrix.h utils.cpp utils.h constants.h FEM_model.cpp FEM_model.h)
