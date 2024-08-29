// Author(s): Matthew Speranza
#ifndef MOLECULAR_DYNAMICS_C_TIMER_H
#define MOLECULAR_DYNAMICS_C_TIMER_H

#include <time.h>

void printClock(clock_t start, clock_t end, char* action);
void printClockDelta(clock_t delta, char* action);
void reportRuntime(clock_t start, clock_t end);

#endif //MOLECULAR_DYNAMICS_C_TIMER_H
