// Author(s): Matthew Speranza
#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H
#include "../system/system.h"

void buildLists(System* system);
void coordToGrid(REAL *xyz, const int nx, const int ny, const int nz);
int gridCoordToGrid(int x, int y, int z, int nx, int ny, int nz);


#endif //NEIGHBORLIST_H
