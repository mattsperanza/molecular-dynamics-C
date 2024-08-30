cmake_minimum_required(VERSION 3.28)
project(molecular_dynamics_C C)
set(CMAKE_C_STANDARD 11)


# Add the executables
add_subdirectory(src/common)
add_subdirectory(src/classical)
add_subdirectory(src/quantum)
include_directories(
        src/classical/include
        src/quantum/include
        src/common/include
)
add_executable(
        molecular_dynamics_C
        ${CLASSICAL}
        ${QUANTUM}
        ${COMMON}
        src/md.c
)

find_library(MATH_LIBRARY m REQUIRED)
if(MATH_LIBRARY)
    target_link_libraries(molecular_dynamics_C PUBLIC ${MATH_LIBRARY})
endif()
find_package(OpenMP REQUIRED)
target_link_libraries(molecular_dynamics_C PRIVATE OpenMP::OpenMP_C)

# Change to O3 to see which loops are vectorized in debug mode
set(FLAGS_DEBUG "-O0;-ffast-math;-fno-math-errno;--verbose;-Wall;--verbose;-fopenmp") # --analyze to run static analysis
set(FLAGS_RELEASE "-O3;-ffast-math;-fno-math-errno;-Wall;-fopenmp")

# Apply compile options to the target
target_compile_options(molecular_dynamics_C PRIVATE "$<$<CONFIG:DEBUG>:${FLAGS_DEBUG}>")
target_compile_options(molecular_dynamics_C PRIVATE "$<$<CONFIG:RELEASE>:${FLAGS_RELEASE}>")