// Author(s): Matthew Speranza
#include <stdio.h>
#include <time.h>

#include "common/include/commandInterpreter.h"
#include "common/include/timer.h"


/**
    * Time the total program execution time and report.
    */
int main(int argc, char* argv[]) {
    clock_t startT, endT;
    startT = clock();
    // Entire mdc program
    commandInterpreter(argc, argv);
    endT = clock();
    reportRuntime(startT, endT);
    return 0;
}

