//
// Created by Matthew Speranza on 8/28/24.
//

#include <stdio.h>
#include "../include/timer.h"

void printClock(double start, double end, char* action){
    double time_sec = end-start;
    long time_min = time_sec / 60;
    long time_hr = time_min / 60;
    long time_day = time_hr / 24;
    if (time_day != 0) {
        printf("Action <%s> runtime day:hr:min:sec --> %-3ld:%-2ld:%-2ld:%-5.3f\n", action, time_day, time_hr, time_min, time_sec);
    } else if (time_hr != 0) {
        printf("Action <%s> runtime hr:min:sec --> %-2ld:%-2ld:%-5.3f\n", action, time_hr, time_min, time_sec);
    } else if (time_min != 0) {
        printf("Action <%s> runtime min:sec --> %-2ld:%-5.3f\n", action, time_min, time_sec);
    } else {
        printf("Action <%s> runtime --> %-8.6lf seconds\n", action, time_sec);
    }
}

void printClockDelta(double dt, char* action){
    double time_sec = dt;
    long time_msec = (long) (time_sec * 1e3);
    long time_min = time_msec / ((long) 1e3 * 60);
    long time_hr = time_min / 60;
    long time_day = time_hr / 24;
    if (time_day != 0) {
        printf("Action <%s> runtime day:hr:min:sec --> %-3ld:%-2ld:%-2ld:%-5.3f\n", action, time_day, time_hr, time_min, time_sec);
    } else if (time_hr != 0) {
        printf("Action <%s> runtime hr:min:sec --> %-2ld:%-2ld:%-5.3f\n", action, time_hr, time_min, time_sec);
    } else if (time_min != 0) {
        printf("Action <%s> runtime min:sec --> %-2ld:%-5.3f\n", action, time_min, time_sec);
    } else {
        printf("Action <%s> runtime --> %-8.6lf seconds\n", action, time_sec);
    }
}

void reportRuntime(double start, double end) {
    double time_sec = end-start;
    long time_min = time_sec / 60;
    long time_hr = time_min / 60;
    long time_day = time_hr / 24;
    if (time_day != 0) {
        printf("\n\nProgram execution time day:hr:min:sec --> %3ld:%2ld:%2ld:%5.3f", time_day, time_hr, time_min, time_sec);
    } else if (time_hr != 0) {
        printf("\n\nProgram execution time hr:min:sec --> %2ld:%2ld:%5.3f", time_hr, time_min, time_sec);
    } else if (time_min != 0) {
        printf("\n\nProgram execution time min:sec --> %2ld:%5.3f", time_min, time_sec);
    } else {
        printf("\n\nProgram execution time --> %-8.6lf seconds", time_sec);
    }
    printf("\n\n");
}
