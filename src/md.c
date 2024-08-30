// Author(s): Matthew Speranza
#include <omp.h>
#include <stdio.h>
#include <time.h>

#include "common/include/commandInterpreter.h"
#include "common/include/timer.h"


/**
    * Time the total program execution time and report.
    */
int main(int argc, char* argv[]) {
    double startT = omp_get_wtime();
    // Entire mdc program
    commandInterpreter(argc, argv);
    double endT = omp_get_wtime();
    reportRuntime(startT, endT);
    return 0;
}

