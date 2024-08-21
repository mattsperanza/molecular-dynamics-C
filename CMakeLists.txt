cmake_minimum_required(VERSION 3.28)
project(molecular_dynamics_C C)
set(CMAKE_C_STANDARD 11)

# Check for Clang existence
find_program(CLANG_PATH clang)
if (NOT CLANG_PATH)
    message(FATAL_ERROR "Clang compiler not found. Please install Clang.")
endif()
set(CMAKE_C_COMPILER ${CLANG_PATH})

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
add_executable(
        commonTest
        src/common/commonTest.c
        ${COMMON}
)

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(molecular_dynamics_C PUBLIC ${MATH_LIBRARY})
endif()

# Change to O3 to see which loops are vectorized in debug mode
set(FLAGS_DEBUG "-O0;-g;-ffast-math;-fno-math-errno;--verbose;-Wall;--verbose;") # --analyze to run static analysis
set(FLAGS_RELEASE "-O3;-ffast-math;-fno-math-errno;-Rpass=loop-vectorize;-Rpass-analysis=loop-vectorize:-Wall;")

# Apply compile options to the target
target_compile_options(molecular_dynamics_C PRIVATE "$<$<CONFIG:DEBUG>:${FLAGS_DEBUG}>")
target_compile_options(molecular_dynamics_C PRIVATE "$<$<CONFIG:RELEASE>:${FLAGS_RELEASE}>")