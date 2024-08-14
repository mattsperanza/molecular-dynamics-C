cmake_minimum_required(VERSION 3.28)
project(molecular_dynamics_C C)

# Set the C standard
set(CMAKE_C_STANDARD 11)

# Check for Clang existence
find_program(CLANG_PATH clang)
if (NOT CLANG_PATH)
    message(FATAL_ERROR "Clang compiler not found. Please install Clang.")
endif()
set(CMAKE_C_COMPILER ${CLANG_PATH})

# Include directories
include_directories(src)
# Add the executable                /main
add_executable(molecular_dynamics_C src/common/script/md.c
        # /classical
        ## /classical/bonded
        ## /classical/forcefields
        ## /classical/nonbonded

        # /common
        ## /common/numerics
        ## /common/parser
        ## /common/script
        src/common/script/commandInterpreter.h src/common/script/commandInterpreter.c
        ## common/utils
        ### common/utils/ds
        src/common/utils/ds/list.h src/common/utils/ds/list.c

        # /quantum
)

# Change to O3 to see which loops are vectorized in debug mode
set(FLAGS_DEBUG "-O0;-ffast-math;-Rpass=loop-vectorize;-Rpass-missed=loop-vectorize;-Rpass-analysis=loop-vectorize")
list(APPEND FLAGS_DEBUG "--verbose;-fno-math-errno")
set(FLAGS_RELEASE "-O3;-ffast-math;-fno-math-errno")

# Apply compile options to the target
target_compile_options(molecular_dynamics_C PRIVATE "$<$<CONFIG:DEBUG>:${FLAGS_DEBUG}>")
target_compile_options(molecular_dynamics_C PRIVATE "$<$<CONFIG:RELEASE>:${FLAGS_RELEASE}>")