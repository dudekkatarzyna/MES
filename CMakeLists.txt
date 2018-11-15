# cmake_minimum_required(VERSION <specify CMake version here>)
project(MES)

set(CMAKE_CXX_STANDARD 14)

add_executable(MES main.cpp Grid/GRID.cpp Grid/GRID.h Grid/Components/NODE.cpp Grid/Components/NODE.h Grid/Components/ELEMENT.cpp Grid/Components/ELEMENT.h Utility/Utility.cpp Utility/Utility.h UniversalElements/UniversalElement.cpp UniversalElements/UniversalElement.h UniversalElements/UniversalElementSides.cpp UniversalElements/UniversalElementSides.h Utility/Calculations.cpp Utility/Calculations.h Utility/Matrix.cpp Utility/Matrix.h UniversalElements/KsiEta.cpp UniversalElements/KsiEta.h)