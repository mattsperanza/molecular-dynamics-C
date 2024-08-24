// Author(s): Matthew Speranza
#include "../include/neighborList.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int comparInt(const void* a, const void* b) {
  return *(int*)a - *(int*)b;
}

// Use bfs to build 1-3, 1-4, and 1-5 pair and path lists
// BFS avoids problems with loops
// Also assign bonded terms
void buildBonded(System* system) {
  printf("Building Bonded Lists\n");
  system->list13 = malloc(sizeof(Vector)*system->nAtoms);
  system->list14 = malloc(sizeof(Vector)*system->nAtoms);
  system->list15 = malloc(sizeof(Vector)*system->nAtoms);
  if(system->list13 == NULL || system->list14 == NULL || system->list15 == NULL) {
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
    visited[i] = calloc(sizeof(int), system->nAtoms);
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
  }
  // 1-5 atoms... this is getting repetitive
  for(int i = 0; i < system->nAtoms; i++) {
    Vector* list15 = vectorCreate(sizeof(int), 10, NULL, INT);
    // Loop over 1-4 atoms
    Vector ilist14 = system->list14[i];
    for(int j = 0; j < ilist14.size; j++) {
      int atomID = ((int*)ilist14.array)[j];
      Vector bfsList = system->list12[atomID];
      for(int k = 0; k < bfsList.size; k++) {
        int atomID2 = ((int*)bfsList.array)[k];
        if(visited[i][atomID2] == 0) {
          vectorAppend(list15, &atomID2); // add the atoms 14 atoms 12 atoms to 15 list
        }
        visited[i][atomID2] = 1;
      }
    }
    free(visited[i]);
  }
  free(visited);

  // Set atom classes
  Vector atoms = system->forceField->atom;
  system->atomClasses = malloc(sizeof(int)*system->nAtoms);
  // Assign atom classes to every atom
  for(int i = 0; i < system->nAtoms; i++) {
    Atom** atomLines = atoms.array;
    int atomType = system->atomTypes[i];
    system->atomClasses[i] = atomLines[atomType-1]->aClass;
  }

  // This was taken from FFX since my way was slightly off
  // Bond paths (loop over bonds and mark atom as part of bond)
  system->bonds = *vectorCreate(sizeof(int*), 300, NULL, INT_PTR);
  system->forceField->bondParams = *vectorCreate(sizeof(Bond*), 300, NULL, OTHER);
  system->piTorsions = *vectorCreate(sizeof(int), 300, NULL, INT);
  system->forceField->piTorsParams = *vectorCreate(sizeof(PiTors*), 300, NULL, OTHER);
  int* partOfBond = calloc(sizeof(int), system->nAtoms); // 0=false
  for(int i = 0; i < system->nAtoms; i++) {
    Vector ilist12 = system->list12[i];
    partOfBond[i] = 1;
    for(int j = 0; j < ilist12.size; j++) {
      int atomID = ((int*)ilist12.array)[j];
      if(partOfBond[atomID] == 1) {
        continue;
      }
      int* bond = malloc(sizeof(int)*2);
      bond[0] = i;
      bond[1] = atomID;
      vectorAppend(&system->bonds, &bond);
      // Search through bonded params for this bond
      Bond** bondParamArray = system->forceField->bond.array;
      bool found = false;
      for(int k = 0; k < system->forceField->bond.size; k++) {
        int* bondClasses = bondParamArray[k]->atomClasses;
        // Get atom classes of bond, sort, and compared to sorted Bond.atomClasses from reading FF file
        int classes[2];
        for(int l = 0; l < 2; l++) {
          classes[l] = system->atomClasses[bond[l]];
        }
        qsort(classes, 2, sizeof(int), comparInt);
        for(int l = 0; l < 2; l++) {
          found = bondClasses[l] == classes[l];
        }
        if(found) {
          vectorAppend(&system->forceField->bondParams, bondParamArray[k]);
          break;
        }
      }
      if(!found) {
        printf("Bonded params not found for bond between classes %d-%d\n", system->atomClasses[bond[0]],
          system->atomClasses[bond[1]]);
        exit(1);
      }
      // TODO: Fix pi-orbital
      // Search Pi-Orbital
      PiTors** piTorsArray = system->forceField->piTors.array;
      found = false;
      for(int k = 0; k < system->forceField->piTors.size; k++) {
        int* piBondClasses = piTorsArray[k]->atomClasses;
        int classes[2];
        for(int l = 0; l < 2; l++) {
          classes[l] = system->atomClasses[bond[l]];
        }
        qsort(classes, 2, sizeof(int), comparInt);
        for(int l = 0; l < 2; l++) {
          found = piBondClasses[l] == classes[l];
        }
        if(found) {
          int index = k;
          vectorAppend(&system->piTorsions, &index);
          vectorAppend(&system->forceField->piTorsParams, piTorsArray[k]);
          break;
        }
      }
    }
  }
  free(partOfBond);

  // Angle paths (loop over bonds and check end bonds)
  system->angles = *vectorCreate(sizeof(int*), 300, NULL, INT_PTR);
  system->forceField->angleParams = *vectorCreate(sizeof(Angle*), 300, NULL, OTHER);
  for(int i = 0; i < system->nAtoms; i++) {
    Vector bonds = system->list12[i];
    int* ids = bonds.array;
    int index = 0;
    for(int j = 0; j < bonds.size; j++) {
      int atomID1 = ids[j];
      index++;
      for(int k = index; k < bonds.size; k++) {
        int atomID2 = ids[k];
        int* angle = malloc(sizeof(int)*3);
        angle[0] = i;
        angle[1] = atomID1;
        angle[2] = atomID2;
        vectorAppend(&system->angles, &angle);
        // Loop through angle params
        Angle** angleParamArray = system->forceField->angle.array;
        int size = system->forceField->angle.size;
        bool found = false;
        for(int l = 0; l < size; l++) {
          int* angleClasses = angleParamArray[l]->aClasses;
          int classes[3];
          for(int m = 0; m < 3; m++) {
            classes[m] = system->atomClasses[angle[m]];
          }
          qsort(classes, 3, sizeof(int), comparInt);
          for(int m = 0; m < 3; m++) {
            found = angleClasses[m] == classes[m];
          }
          if(found) {
            vectorAppend(&system->forceField->angleParams, angleParamArray[l]);
            break;
          }
        }
        if(!found) {
          printf("Failed to find angle params for angle between classes %d-%d-%d\n", system->atomClasses[angle[0]],
            system->atomClasses[angle[1]], system->atomClasses[angle[2]]);
          exit(1);
        }
      }
    }
  }

  // StrBend & UrayBrad
  Vector angles = system->angles;
  system->strBendAngleIndex = *vectorCreate(sizeof(int), 300, NULL, INT);
  system->forceField->strBendParams = *vectorCreate(sizeof(StrBend*), 300, NULL, OTHER);
  system->urayBradIndex = *vectorCreate(sizeof(int), 300, NULL, INT);
  system->forceField->urayBradParams = *vectorCreate(sizeof(UReyBrad*), 300, NULL, OTHER);
  for(int i = 0; i < angles.size; i++) {
    int* angle = ((int**)angles.array)[i];
    // Check for stretch bend
    // TODO: Fix the issue with sorting
    StrBend** strbendParamArray = system->forceField->strBend.array;
    int size = system->forceField->strBend.size;
    bool found = false;
    for(int j = 0; j < size; j++) {
      int* strBendClasses = strbendParamArray[j]->atomClasses;
      int classes[3];
      for(int k = 0; k < 3; k++) {
        classes[k] = system->atomClasses[angle[k]];
      }
      qsort(classes, 3, sizeof(int), comparInt);
      for(int k = 0; k < 3; k++) {
        found = strBendClasses[k] == classes[k];
      }
      if(found) {
        int index = i;
        vectorAppend(&system->strBendAngleIndex, &index);
        vectorAppend(&system->forceField->strBendParams, strbendParamArray[j]);
        break;
      }
    }
    // Check for Uray-Bradley
    UReyBrad** urayBradParamArray = system->forceField->uRayBrad.array;
    size = system->forceField->uRayBrad.size;
    found = false;
    for(int j = 0; j < size; j++) {
      int* urayBradClasses = urayBradParamArray[j]->atomClasses;
      int classes[3];
      for(int k = 0; k < 3; k++) {
        classes[k] = system->atomClasses[angle[k]];
      }
      qsort(classes, 3, sizeof(int), comparInt);
      for(int k = 0; k < 3; k++) {
        found = urayBradClasses[k] == classes[k];
      }
      if(found) {
        int index = i;
        vectorAppend(&system->urayBradIndex, &index);
        vectorAppend(&system->forceField->urayBradParams, &index);
        break;
      }
    }
  }

  // Torsion paths (loop over angles and check end bonds)
  system->torsions = *vectorCreate(sizeof(int*), 300, NULL, INT_PTR);
  system->forceField->torsionParams = *vectorCreate(sizeof(Torsion*), 300, NULL, OTHER);
  Vector bonds = system->bonds;
  for(int i = 0; i < bonds.size; i++) {
    int* bond = ((int**)bonds.array)[i];
    int atomID1 = bond[0];
    int atomID2 = bond[1];
    Vector atom1Bonds = system->list12[atomID1];
    Vector atom2Bonds = system->list12[atomID2];
    for(int j = 0; j < atom1Bonds.size; j++) {
      int atomID3 = ((int*)atom1Bonds.array)[j];
      if(atomID3 != atomID1 && atomID3 != atomID2) {
        for(int k = 0; k < atom2Bonds.size; k++) {
          int atomID4 = ((int*)atom2Bonds.array)[k];
          if(atomID4 != atomID1 && atomID4 != atomID2) {
            int* torsion = malloc(sizeof(int)*4);
            torsion[0] = atomID1;
            torsion[1] = atomID2;
            torsion[2] = atomID3;
            torsion[3] = atomID4;
            vectorAppend(&system->torsions, &torsion);
            // Match torsion params
            Torsion** torsionParamArray = system->forceField->torsion.array;
            int size = system->forceField->torsion.size;
            bool found = false;
            for(int l = 0; l < size; l++) {
              int* torsionClasses = torsionParamArray[l]->atomClasses;
              int classes[4];
              for(int m = 0; m < 4; m++) {
                classes[m] = system->atomClasses[torsion[m]];
              }
              qsort(classes, 4, sizeof(int), comparInt);
              for(int m = 0; m < 4; m++) {
                found = torsionClasses[m] == classes[m];
              }
              if(found) {
                vectorAppend(&system->forceField->torsionParams, torsionParamArray[l]);
                break;
              }
            }
            if(!found) {
              printf("Failed to find torsion parameters for class set %d-%d-%d-%d\n", system->atomClasses[torsion[0]],
                system->atomClasses[torsion[1]], system->atomClasses[torsion[2]], system->atomClasses[torsion[3]]);
              exit(1);
            }
          }
        }
      }
    }
  }

  // Torsion-Torsion - TODO: Fix this to match ffx
  system->torsionTorsion = *vectorCreate(sizeof(int*), 300, NULL, INT_PTR);
  for(int i = 0; i < angles.size; i++) {
    int* angle = ((int**)angles.array)[i];
    int atomID1 = angle[0];
    int atomID2 = angle[1];
    int atomID3 = angle[2];
    Vector atom1Bonds = system->list12[atomID1];
    Vector atom3Bonds = system->list12[atomID3];
    for(int j = 0; j < atom1Bonds.size; j++) {
      int atomID0 = ((int*)atom1Bonds.array)[j];
      if(atomID0 != atomID2 && atomID0 != atomID3) {
        for(int k = 0; k < atom3Bonds.size; k++) {
          int atomID4 = ((int*)atom3Bonds.array)[k];
          if(atomID4 != atomID0 && atomID4 != atomID1 && atomID4 != atomID2) {
            int* tortors = malloc(sizeof(int)*5);
            tortors[0] = atomID0;
            tortors[1] = atomID1;
            tortors[2] = atomID2;
            tortors[3] = atomID3;
            tortors[4] = atomID4;
            vectorAppend(&system->torsionTorsion, &tortors);
            // TODO: Match torsion torsion type
          }
        }
     }
    }
  }
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
  REAL rCut2 = system->vdwCutoff + system->realspaceBuffer;
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
  REAL rCut2 = system->vdwCutoff + system->realspaceBuffer;
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
  REAL rCut = system->vdwCutoff + system->realspaceBuffer;
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
