// Author(s): Matthew Speranza
#include "../include/neighborList.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Use bfs to build 1-3 and 1-4 lists
void buildBonded(System* system) {
  system->list13 = malloc(sizeof(Vector)*system->nAtoms);
  system->list14 = malloc(sizeof(Vector)*system->nAtoms);
  if(system->list13 == NULL || system->list14 == NULL) {
    printf("Failed to allocate memory for list13 or list14 in buildBonded\n");
    exit(1);
  }
  int** visited = malloc(sizeof(int*)*system->nAtoms);
  if(visited == NULL) {
    printf("Failed to allocate memory for visited array in buildBonded\n");
    exit(1);
  }
  // 1-3 atoms
  for(int i = 0; i < system->nAtoms; i++) {
    Vector* list13 = vectorCreate(sizeof(int), 10, NULL, INT);
    // Mark all atoms as unvisited - 0=false
    visited[i] = malloc(sizeof(int)*system->nAtoms);
    if(visited[i] == NULL) {
      printf("Failed to allocate memory for visited array in buildBonded\n");
      exit(1);
    }
    // Mark unwanted atoms as visited
    visited[i][i] = 1; // the atom itself
    Vector ilist12 = system->list12[i];
    for(int j = 0; j < ilist12.size; j++) {
      int atomID = ((int*)ilist12.array)[j];
      visited[i][atomID] = 1; // the atoms 12 atoms
      Vector bfsList = system->list12[atomID];
      for(int k = 0; k < bfsList.size; k++) {
        int atomID2 = ((int*)bfsList.array)[k];
        if(visited[i][atomID2] == 0) {
          vectorAppend(list13, &atomID2); // add the atoms 12 atoms 12 atoms to 13 list
        }
        visited[i][atomID2] = 1; // the atoms 13 atoms
      }
    }
    system->list13[i] = *list13;
  }
  // 1-4 atoms
  for(int i = 0; i < system->nAtoms; i++) {
    Vector* list14 = vectorCreate(sizeof(int), 10, NULL, INT);
    // Loop over 1-3 atoms
    Vector ilist13 = system->list13[i];
    for(int j = 0; j < ilist13.size; j++) {
      int atomID = ((int*)ilist13.array)[j];
      Vector bfsList = system->list12[atomID];
      for(int k = 0; k < bfsList.size; k++) {
        int atomID2 = ((int*)bfsList.array)[k];
        if(visited[i][atomID2] == 0) {
          vectorAppend(list14, &atomID2); // add the atoms 13 atoms 12 atoms to 14 list
        }
        visited[i][atomID2] = 1; // mark the atoms 14 atoms
      }
    }
    system->list14[i] = *list14;
    free(visited[i]);
  }
  free(visited);
}

int indexGrid(int x, int y, int z, int nx, int ny, int nz) {
  // Shift the index to the correct cell inside box
  if(x >= nx) { x -= nx; }
  else if(x < 0) { x += nx; }
  if(y >= ny) { y -= ny; }
  else if(y < 0) { y += ny; }
  if(z >= nz) { z -= nz; }
  else if(z < 0) { z += nz; }
  int index = nx * x + ny * y + z;
  assert(index >= 0 && index <= nx*ny*nz);
  return index;
}

/**
 * Applies the minimum image convention to a distance vector
 * @return a new distance
 */
extern REAL imageDx(REAL dx, REAL axisLen) {
 while(dx > axisLen/2 || dx <= -axisLen/2)
  dx = dx > 0 ? dx - axisLen : dx + axisLen;
 return dx;
};

void addCellToList(Vector* cell, Vector* list, System* system, int atomID, REAL aLen, REAL bLen, REAL cLen) {
  REAL* pos = system->X;
  REAL rCut2 = system->realspaceCutoff + system->realspaceBuffer;
  rCut2 *= rCut2;
  for(int i = 0; i < cell->size; i++) {
    int atomID2 = ((int*)cell->array)[i];
    REAL dx = imageDx(pos[atomID*3] - pos[atomID2*3], aLen);
    REAL dy = imageDx(pos[atomID*3+1] - pos[atomID2*3+1], bLen);
    REAL dz = imageDx(pos[atomID*3+2] - pos[atomID2*3+2], cLen);
    REAL r2 = dx*dx + dy*dy + dz*dz;
    if(r2 < rCut2) {
      vectorAppend(list, &atomID2);
    }
  }
}

void addCellToListSelf(Vector* cell, Vector* list, System* system, int atomID, REAL aLen, REAL bLen, REAL cLen) {
  REAL* pos = system->X;
  REAL rCut2 = system->realspaceCutoff + system->realspaceBuffer;
  rCut2 *= rCut2;
  for(int i = 0; i < cell->size; i++) {
    int atomID2 = ((int*)cell->array)[i];
    REAL dx = imageDx(pos[atomID*3] - pos[atomID2*3], aLen);
    REAL dy = imageDx(pos[atomID*3+1] - pos[atomID2*3+1], bLen);
    REAL dz = imageDx(pos[atomID*3+2] - pos[atomID2*3+2], cLen);
    REAL r2 = dx*dx + dy*dy + dz*dz;
    if(r2 < rCut2 && atomID <= atomID2) {
      vectorAppend(list, &atomID2);
    }
  }
}

void buildVerlet(System* system) {
  // Find axis lengths and grid spacing
  REAL* a = system->boxDim[0];
  // len(vec) = sqrt(dot(vec, vec))
  REAL aLen = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  REAL* b = system->boxDim[1];
  REAL bLen = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
  REAL* c = system->boxDim[2];
  REAL cLen = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  // a dot (b cross c) = volume
  system->volume = a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0]) + a[2]*(b[0]*c[1] - b[1]*c[0]);
  system->particleDensity = system->nAtoms / system->volume;
  system->verletList = malloc(sizeof(Vector)*system->nAtoms);
  system->realspaceBuffer = 2;
  float num = 16 / (aLen + bLen + cLen);
  // Set number of grid cells in each direction
  int nX = num * aLen + 1;
  if(nX <= 1) {
    nX = 2;
  }
  float xCubeLen = aLen / nX;
  int nY = num * bLen + 1;
  if(nY <= 1) {
    nY = 2;
  }
  float yCubeLen = bLen / nY;
  int nZ = num * cLen + 1;
  if(nZ <= 1) {
    nZ = 2;
  }
  float zCubeLen = cLen / nZ;
  int nCells = nX * nY * nZ;
  // Loop over all atoms and assign them to grid cells
  Vector* grid = calloc(sizeof(Vector), nX*nY*nZ);
  for(int i = 0; i < system->nAtoms; i++) {
    REAL x = system->X[i*3] - system->minDim[0]; // shift unit cell into +x, +y, +z octant
    REAL y = system->X[i*3+1] - system->minDim[1];
    REAL z = system->X[i*3+2] - system->minDim[2];
    int gridX = x / xCubeLen;
    int gridY = y / yCubeLen;
    int gridZ = z / zCubeLen;
    int index = indexGrid(gridX, gridY, gridZ, nX, nY, nZ);
    if(grid[index].array == NULL) {
      grid[index] = *vectorCreate(sizeof(int), 32, NULL, INT);
    }
    int atomID = i;
    vectorAppend(&grid[index], &atomID);
  }
  // Loop over all atoms again and loop over half of the neighboring cells to build list
  REAL rCut = system->realspaceCutoff + system->realspaceBuffer;
  long interactionsCell = 0;
  int* visitedCells = calloc(sizeof(int), nCells);
  for(int i = 0; i < system->nAtoms; i++) {
    system->verletList[i] = *vectorCreate(sizeof(int), 1e3, NULL, INT);
    REAL x = system->X[i*3] - system->minDim[0];
    REAL y = system->X[i*3+1] - system->minDim[1];
    REAL z = system->X[i*3+2] - system->minDim[2];
    int gridX = x / xCubeLen;
    int searchX = rCut/xCubeLen+1;
    int gridY = y / yCubeLen;
    int searchY = rCut/yCubeLen+1;
    int gridZ = z / zCubeLen;
    int searchZ = rCut/zCubeLen+1;
    // Add atoms from the current cell
    int cellID = indexGrid(gridX, gridY, gridZ, nX, nY, nZ);
    visitedCells[cellID] = 1;
    addCellToListSelf(&grid[cellID], &system->verletList[i], system, i, aLen, bLen, cLen);
    // Add atoms from y-direction line of cells
    for(int j = gridY+1; j <= gridY+searchY; j++) {
      int index = indexGrid(gridX, j, gridZ, nX, nY, nZ);
      if(visitedCells[index] != 1) {
        addCellToList(&grid[index], &system->verletList[i], system, i, aLen, bLen, cLen);
      }
      visitedCells[index] = 1;
    }
    // Add atoms from z-direction half-plane of cells
    for(int j = gridY-searchY; j <= gridY+searchY; j++) {
      for(int k = gridZ+1; k <= gridZ+searchZ; k++) {
        int index = indexGrid(gridX, j, k, nX, nY, nZ);
        if(visitedCells[index] != 1) {
          addCellToList(&grid[index], &system->verletList[i], system, i, aLen, bLen, cLen);
        }
        visitedCells[index] = 1;
      }
    }
    // Add atoms from x-direction half-cube of cells
    for(int j = gridY-searchY; j <= gridY+searchY; j++) {
      for(int k = gridZ-searchZ; k <= gridZ+searchZ; k++) {
        for(int l = gridX+1; l <= gridX+searchX; l++) {
          int index = indexGrid(j, k, l, nX, nY, nZ);
          if(visitedCells[index] != 1) {
            addCellToList(&grid[index], &system->verletList[i], system, i, aLen, bLen, cLen);
          }
          visitedCells[index] = 1;
        }
      }
    }
    interactionsCell+=system->verletList[i].size;
    memset(visitedCells, 0, sizeof(int)*nCells);
  }
  free(visitedCells);
  //printf("Interactions: %ld\n", interactionsCell);
  //printf("Vector: ");
  //vectorPrint(&system->verletList[0], system->verletList[0].size);
  for(int i = 0; i < nCells; i++) {
    vectorBackingFree(&grid[i]);
  }
  free(grid);

}


void buildLists(System* system) {
  buildBonded(system);
  buildVerlet(system);
};
