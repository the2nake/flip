cmake_minimum_required(VERSION 3.10)

project(
  FlipFluidSimulator
  LANGUAGES C)

find_package(SDL3 REQUIRED CONFIG REQUIRED COMPONENTS SDL3-shared)

add_link_options("-lm")
add_executable(ffs src/main.c)

target_link_libraries(ffs PRIVATE SDL3::SDL3)
set_target_properties(ffs PROPERTIES
  C_STANDARD 23
  C_STANDARD_REQUIRED TRUE
)

# target_compile_options(ffs PRIVATE -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl)
