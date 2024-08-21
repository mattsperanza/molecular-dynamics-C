// Author(s): Matthew Speranza
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector.h>
#include "../include/xyz.h"

#include <limits.h>

int splitLine(char* line, char delim, char**);

/**
 * XYZ Format (Tinker/AMOEBA):
 * {nAtoms} {remark}
 * 1 AtomString x y z atomType bondedAtomID bondedAtomID ...
 * 2 AtomString x y z atomType bondedAtomID bondedAtomID ...
 * ...
 * nAtoms AtomString x y z atomType bondedAtomID bondedAtomID ...
 *
 * @param system system to fill out
 */
void readXYZ(System* system, char* structureFileName) {
 FILE* f = fopen(structureFileName, "r");
 if(f == NULL) {
  printf("Failed to open file %s", structureFileName);
  exit(1);
 }
 system->structureFileName = structureFileName;
 system->patchFiles = *vectorCreate(sizeof(char*), 1, NULL, CHAR_PTR);
 system->forceFieldFile = malloc(sizeof(char)*1000);
 // Read whitespace and first line
 int lineSize = 1e3;
 char line[lineSize];
 int nAtoms = -1;
 if(fgets(line, lineSize, f) == NULL) {
  printf("Failed to read from file: %s", structureFileName);
  exit(1);
 }
 while(strcmp(&line[0], "\n") == 0) {
  if(fgets(line, lineSize, f) == NULL) {
   printf("Failed to read from file: %s", structureFileName);
   exit(1);
  }
 }
 if(*line == EOF) {
  printf("Error reading line from file: %s", structureFileName);
  exit(1);
 }
 system->remark = strdup(line);
 nAtoms = atoi(strtok(line, " "));
 if(nAtoms == -1) {
  printf("Failed to find number of atoms in file %s", structureFileName);
  exit(1);
 }
 system->nAtoms = nAtoms;
 // 2d arrays
 system->multipoles = malloc(sizeof(REAL*)*nAtoms);
 system->list12 = malloc(sizeof(Vector)*nAtoms);
 system->atomNames = malloc(sizeof(char*)*nAtoms);
 // 1d arrays
 system->atomTypes = malloc(sizeof(int)*nAtoms*3);
 system->X = malloc(sizeof(REAL)*nAtoms*3);
 system->M = malloc(sizeof(REAL)*nAtoms);
 system->V = malloc(sizeof(REAL)*nAtoms*3);
 system->A = malloc(sizeof(REAL)*nAtoms*3);
 system->F = malloc(sizeof(REAL)*nAtoms*3);
 system->lambdas = malloc(sizeof(REAL)*nAtoms*3);
 system->protons = malloc(sizeof(REAL)*nAtoms*3);
 system->valence = malloc(sizeof(REAL)*nAtoms*3);
 // Spatial
 for(int i = 0; i < 3; i++) {
  system->minDim[i] = INT_MAX;
  for(int j = 0; j < 3; j++) {
   system->boxDim[i][j] = -1.0f; // check later to see if box dim was set
  }
 }
 system->pmeGridspace = malloc(sizeof(int)*3);
 // Read atom lines
 Vector* bonded = vectorCreate(sizeof(int), 10, NULL, INT);
 for(int i = 0; i < nAtoms; i++) {
  if(fgets(line, lineSize, f) == NULL) {
   printf("Failed to read on line %d of %s!", i, structureFileName);
   exit(1);
  }
  int atomIndex = atoi(strtok(line, " "))-1;
  system->lambdas[atomIndex] = 1.0;
  system->atomNames[atomIndex] = strdup(strtok(NULL, " "));
  system->X[atomIndex*3] = atof(strtok(NULL, " "));
  if(system->X[atomIndex*3] < system->minDim[0]) {
   system->minDim[0] = system->X[atomIndex*3];
  }
  system->X[atomIndex*3+1] = atof(strtok(NULL, " "));
  if(system->X[atomIndex*3+1] < system->minDim[1]) {
   system->minDim[1] = system->X[atomIndex*3+1];
  }
  system->X[atomIndex*3+2] = atof(strtok(NULL, " "));
  if(system->X[atomIndex*3+2] < system->minDim[2]) {
   system->minDim[2] = system->X[atomIndex*3+2];
  }
  system->atomTypes[atomIndex] = atoi(strtok(NULL, " "));
  char* str = strtok(NULL, " ");
  while(str != NULL) {
   int bondedID = atoi(str)-1; // allocate new address
   vectorAppend(bonded, &bondedID);
   str = strtok(NULL, " ");
  }
  system->list12[atomIndex] = *vectorCopy(bonded);
  vectorClear(bonded);
  assert(atomIndex == i);
 }
 //printXYZ(system);
 vectorBackingFree(bonded);
 fclose(f);
};

void printXYZ(System* system) {
 assert(system != NULL);
 assert(system->X != NULL);
 for(int i = 0; i < system->nAtoms; i++) {
  double x = system->X[i*3];
  double y = system->X[i*3+1];
  double z = system->X[i*3+2];
  printf("Atom %d Name %s Type %d R=(%lf,%lf,%lf) Bonded=[", i+1, system->atomNames[i], system->atomTypes[i], x, y, z);
  Vector bonded = system->list12[i];
  for(int j = 0; j < bonded.size; j++) {
   int bondedAtomID = ((int*)bonded.array)[j];
   printf("%d,", bondedAtomID);
  }
  printf("]\n\n");
 }
}

void readARC(System* system, char* structureName);
void writeXYZFile(System* system, char* outputFileName);
