// Author(s): Matthew Speranza
#ifndef XYZ_H
#define XYZ_H
#include "../system/system.h"

/**
 * This provides functions to read and write xyz files and fill out/push out the system state as fully as possible.
 */
void readXYZ(System* system, char* structureFileName);
void readARC(System* system, char* structureFileName);
void printXYZ(System* system);
void writeXYZ(System* system, char* outputFile);

#endif //XYZ_H
