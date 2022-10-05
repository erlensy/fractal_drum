cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(fractal_drum)

include_directories(src/fractal.h src/lattice.h)

find_package(Armadillo REQUIRED)

include_directories(${ARMADILLO_INCLUDE_DIRS})

file(GLOB SOURCES fractal_drum src/*.cpp src/*.h src/*.py)

add_executable(fractal_drum ${SOURCES})

target_link_libraries(fractal_drum ${ARMADILLO_LIBRARIES})
