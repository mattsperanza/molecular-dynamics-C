set(PWD src/common/)
set(
        COMMON
        ${PWD}system/defines.h # First file included
        ${PWD}system/system.h
        # Paths from root
        # numerics/
        ${PWD}/numerics/neighborList.c
        # parsers/
        ${PWD}parsers/forceFieldReader.c
        ${PWD}parsers/keyReader.c
        ${PWD}parsers/xyz.c
        # scripts/
        ${PWD}scripts/commandInterpreter.c
        # system/
        ${PWD}system/system.h
        # utils/
        ## utils/ds
        ${PWD}utils/ds/vector.c
        PARENT_SCOPE
)
