// Author(s): Matthew Speranza
#ifndef MOLECULAR_DYNAMICS_C_TIMER_H
#define MOLECULAR_DYNAMICS_C_TIMER_H

#include <time.h>

void printClock(double start, double end, char* action);
void printClockDelta(double delta, char* action);
void reportRuntime(double start, double end);

#endif //MOLECULAR_DYNAMICS_C_TIMER_H
