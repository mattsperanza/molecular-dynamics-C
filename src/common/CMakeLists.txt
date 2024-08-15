set(CURRENT_DIR src/common/)
set(
        COMMON
        # Paths from root
        # numerics/
        # parsers/
        # script/
        ${CURRENT_DIR}scripts/commandInterpreter.c
        # utils/
        ## utils/ds
        ${CURRENT_DIR}utils/ds/linkedList.c
        ${CURRENT_DIR}utils/ds/vector.c
        PARENT_SCOPE
)
include_directories(include)
add_executable(
        commonTest
        tests/commonTest.c
)
