set(PWD src/classical/)
set(
        CLASSICAL
        # bonded/
        ${PWD}bonded/bonded.c
        # forcefields/
        # nonbonded/
        ${PWD}nonbonded/direct.c
        PARENT_SCOPE
)