// Author(s): Matthew Speranza
#include "../include/neighborList.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Use bfs to build 1-3 and 1-4 lists
void buildBonded(System* system) {
  printf("Building Bonded Lists\n");
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

/**
 * Moves xGridCoord into xyz and returns the index of the grid cell.
 * @param xyz
 * @param nx
 * @param ny
 * @param nz
 * @return
 */
int coordToGrid(REAL* xyz, const int nx, const int ny, const int nz) {
  // Shift inside box
  REAL x = xyz[0];
  REAL y = xyz[1];
  REAL z = xyz[2];
  while (x >= nx || x < 0) { x = x<0 ? x + 1.0 : x - 1.0; }
  while (y >= ny || y < 0) { y = y<0 ? y + 1.0 : y - 1.0; }
  while (z >= nz || z < 0) { z = z<0 ? z + 1.0 : z - 1.0; }
  int xCoord = nx * x; // Round down
  int yCoord = ny * y;
  int zCoord = nz * z;
  int index = gridCoordToGrid(xCoord, yCoord, zCoord, nx, ny, nz);
  xyz[0] = xCoord;
  xyz[1] = yCoord;
  xyz[2] = zCoord;
  assert(index >= 0 && index <= nx*ny*nz);
  return index;
}

int gridCoordToGrid(int x, int y, int z, const int nx, const int ny, const int nz) {
  if(x < 0 || x >= nx) {
    x = x<0 ? x + nx : x - nx;
  }
  if(y < 0 || y >= ny) {
    y = y<0 ? y + ny : y - ny;
  }
  if (z < 0 || z >= nz) {
    z = z<0 ? z + nz : z - nz;
  }
  // Unique index for each grid cell
  return z + nz * (y + ny * x);
}

/**
 * Applies the minimum image convention to a distance vector
 * @return a new distance
 */
extern REAL imageDx(REAL dx, const REAL axisLen) {
 while(dx > axisLen/2 || dx <= -axisLen/2)
  dx = dx > 0 ? dx - axisLen : dx + axisLen;
 return dx;
};

void addCellToList(const Vector* cell, Vector* list, const System* system,
  const int atomID, const REAL aLen, const REAL bLen, const REAL cLen) {
  REAL* pos = system->X;
  REAL rCut2 = system->realspaceCutoff + system->realspaceBuffer;
  rCut2 *= rCut2;
  for(int i = 0; i < cell->size; i++) {
    int atomID2 = ((int*)cell->array)[i];
    REAL dx = imageDx(pos[atomID*3] - pos[atomID2*3], aLen);
    REAL dy = imageDx(pos[atomID*3+1] - pos[atomID2*3+1], bLen);
    REAL dz = imageDx(pos[atomID*3+2] - pos[atomID2*3+2], cLen);
    REAL r2 = dx*dx + dy*dy + dz*dz;
    if(r2 <= rCut2) {
      vectorAppend(list, &atomID2);
    }
  }
}

void addCellToListSelf(const Vector* cell, Vector* list, const System* system, const int atomID,
  const REAL aLen, const REAL bLen, const REAL cLen) {
  REAL* pos = system->X;
  REAL rCut2 = system->realspaceCutoff + system->realspaceBuffer;
  rCut2 *= rCut2;
  for(int i = 0; i < cell->size; i++) {
    int atomID2 = ((int*)cell->array)[i];
    REAL dx = imageDx(pos[atomID*3] - pos[atomID2*3], aLen);
    REAL dy = imageDx(pos[atomID*3+1] - pos[atomID2*3+1], bLen);
    REAL dz = imageDx(pos[atomID*3+2] - pos[atomID2*3+2], cLen);
    REAL r2 = dx*dx + dy*dy + dz*dz;
    if(r2 <= rCut2 && atomID < atomID2) {
      vectorAppend(list, &atomID2);
    }
  }
}

void toFractional(REAL* xyz, REAL boxDim[3][3]) {
  REAL x = xyz[0];
  REAL y = xyz[1];
  REAL z = xyz[2];
  xyz[0] = x*boxDim[0][0] + y*boxDim[1][0] + z*boxDim[2][0];
  xyz[1] = x*boxDim[0][1] + y*boxDim[1][1] + z*boxDim[2][1];
  xyz[2] = x*boxDim[0][2] + y*boxDim[1][2] + z*boxDim[2][2];
}

void invert3x3(REAL boxDim[3][3], REAL boxDimInv[3][3]) {
  REAL det = boxDim[0][0]*(boxDim[1][1]*boxDim[2][2] - boxDim[1][2]*boxDim[2][1]) -
    boxDim[0][1]*(boxDim[1][0]*boxDim[2][2] - boxDim[1][2]*boxDim[2][0]) +
      boxDim[0][2]*(boxDim[1][0]*boxDim[2][1] - boxDim[1][1]*boxDim[2][0]);
  boxDimInv[0][0] = (boxDim[1][1]*boxDim[2][2] - boxDim[1][2]*boxDim[2][1]) / det;
  boxDimInv[0][1] = (boxDim[0][2]*boxDim[2][1] - boxDim[0][1]*boxDim[2][2]) / det;
  boxDimInv[0][2] = (boxDim[0][1]*boxDim[1][2] - boxDim[0][2]*boxDim[1][1]) / det;
  boxDimInv[1][0] = (boxDim[1][2]*boxDim[2][0] - boxDim[1][0]*boxDim[2][2]) / det;
  boxDimInv[1][1] = (boxDim[0][0]*boxDim[2][2] - boxDim[0][2]*boxDim[2][0]) / det;
  boxDimInv[1][2] = (boxDim[0][2]*boxDim[1][0] - boxDim[0][0]*boxDim[1][2]) / det;
  boxDimInv[2][0] = (boxDim[1][0]*boxDim[2][1] - boxDim[1][1]*boxDim[2][0]) / det;
  boxDimInv[2][1] = (boxDim[0][1]*boxDim[2][0] - boxDim[0][0]*boxDim[2][1]) / det;
  boxDimInv[2][2] = (boxDim[0][0]*boxDim[1][1] - boxDim[0][1]*boxDim[1][0]) / det;
}

void buildVerlet(System* system) {
  printf("Building Neighbor List\n");
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
  float num = 64 / (aLen + bLen + cLen);
  // Set number of grid cells in each direction
  int nX = num * aLen + 1; // round up
  if(nX <= 1) { // minimum of 2 cells
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
  // Calculate the inverse of the box dim matrix
  REAL boxDimInv[3][3];
  invert3x3(system->boxDim, boxDimInv);
  // Loop over all atoms and assign them to grid cells
  Vector* grid = calloc(sizeof(Vector), nX*nY*nZ);
  for(int i = 0; i < system->nAtoms; i++) {
    REAL xyz[3] = {system->X[i*3], system->X[i*3+1], system->X[i*3+2]};
    toFractional(xyz, boxDimInv);
    int index = coordToGrid(xyz, nX, nY, nZ);
    if(grid[index].array == NULL) {
      grid[index] = *vectorCreate(sizeof(int), 32, NULL, INT);
    }
    int atomID = i;
    vectorAppend(&grid[index], &atomID);
  }
  // Loop over all atoms again and loop over half of the neighboring cells to build list
  REAL rCut = system->realspaceCutoff + system->realspaceBuffer;
  int searchX = rCut / xCubeLen + 1; // round up
  int searchY = rCut / yCubeLen + 1;
  int searchZ = rCut / zCubeLen + 1;
  long interactionsCell = 0;
  for(int i = 0; i < system->nAtoms; i++) {
    system->verletList[i] = *vectorCreate(sizeof(int), 1e3, NULL, INT);
    REAL xyz[3] = {system->X[i*3], system->X[i*3+1], system->X[i*3+2]};
    toFractional(xyz, boxDimInv);
    int cellID = coordToGrid(xyz, nX, nY, nZ);
    int gridX = xyz[0];
    int gridY = xyz[1];
    int gridZ = xyz[2];
    addCellToListSelf(&grid[cellID], &system->verletList[i], system, i, aLen, bLen, cLen);
    // z-direction loops are best for cache?
    // Add atoms from z-direction line of cells (excluding self)
    for(int j = gridZ+1; j <= gridZ+searchZ; j++) { // Half z-range
      int index = gridCoordToGrid(gridX, gridY, j, nX, nY, nZ);
      addCellToList(&grid[index], &system->verletList[i], system, i, aLen, bLen, cLen);
    }
    // Add atoms from yz-direction half-plane of cells
    for(int j = gridY+1; j <= gridY+searchY; j++) { // Half y-range
      for(int k = gridZ-searchZ; k <= gridZ+searchZ; k++) { // Full z-range
        int index = gridCoordToGrid(gridX, j, k, nX, nY, nZ);
        addCellToList(&grid[index], &system->verletList[i], system, i, aLen, bLen, cLen);
      }
    }
    // Add atoms from x-direction half-cube of cells
    for(int j = gridX+1; j <= gridX+searchX; j++) { // Half x-range
      for(int k = gridY-searchY; k <= gridY+searchY; k++) { // Full y-range
        for(int l = gridZ-searchZ; l <= gridZ+searchZ; l++) { // Full z-range
          int index = gridCoordToGrid(j, k, l, nX, nY, nZ);
          addCellToList(&grid[index], &system->verletList[i], system, i, aLen, bLen, cLen);
        }
      }
    }
    interactionsCell+=system->verletList[i].size;
  }
  printf("Interactions: %ld\n", interactionsCell);
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
