# cmake_minimum_required(VERSION <specify CMake version here>)
project(MES)

set(CMAKE_CXX_STANDARD 14)

add_executable(MES main.cpp GRID.cpp GRID.h NODE.cpp NODE.h ELEMENT.cpp ELEMENT.h Utility.cpp Utility.h UniversalElement.cpp UniversalElement.h)