// Author(s): Matthew Speranza
#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H
#include "../system/system.h"

void buildLists(System* system);
int coordToGrid(REAL* xyz, int nx, int ny, int nz);
int gridCoordToGrid(int x, int y, int z, int nx, int ny, int nz);


#endif //NEIGHBORLIST_H
