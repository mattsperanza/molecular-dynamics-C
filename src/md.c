// Author(s): Matthew Speranza
#include <stdio.h>
#include <time.h>

#include "common/include/commandInterpreter.h"

void reportRuntime(struct timespec start, struct timespec end) {
    long time_nsec = end.tv_nsec - start.tv_nsec;
    long time_msec = time_nsec / (long) 1e6;
    float time_sec = (double) time_msec / 1e3;
    long time_min = time_msec / ((long) 1e3 * 60);
    long time_hr = time_min / 60;
    long time_day = time_hr / 24;
    if (time_day != 0) {
        printf("\n\nExecution time day:hr:min:sec --> %3ld:%2ld:%2ld:%5.3f", time_day, time_hr, time_min, time_sec);
    } else if (time_hr != 0) {
        printf("\n\nExecution time hr:min:sec --> %2ld:%2ld:%5.3f", time_hr, time_min, time_sec);
    } else if (time_min != 0) {
        printf("\n\nExecution time min:sec --> %2ld:%5.3f", time_min, time_sec);
    } else {
        printf("\n\nExecution time seconds --> %7.5f", time_sec);
    }
    printf("\n\n");
}

/**
    * Time the total program execution time and report.
    */
int main(int argc, char* argv[]) {
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    commandInterpreter(argc, argv);
    clock_gettime(CLOCK_MONOTONIC, &end);
    reportRuntime(start, end);
    return 0;
}

