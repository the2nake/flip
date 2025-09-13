cmake_minimum_required(VERSION 3.10)

project(
  FlipFluidSimulator
  LANGUAGES CXX)

find_package(SDL3 REQUIRED CONFIG REQUIRED COMPONENTS SDL3-shared)

add_executable(ffs src/main.cpp)

target_link_libraries(ffs PRIVATE SDL3::SDL3)
