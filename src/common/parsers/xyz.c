#include "../include/xyz.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector.h>

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
  printf("Couldn't open file %s", structureFileName);
  exit(1);
 }
 system->structureFileName = structureFileName;
 // Read whitespace and first line
 int lineSize = 1e3;
 char* line = malloc(sizeof(char)*lineSize);
 int nAtoms = -1;
 fgets(line, lineSize, f);
 while(strcmp(&line[0], "\n") == 0) {
  fgets(line, lineSize, f);
  if(line == NULL) {
   printf("Failed to read from file: %s", structureFileName);
   exit(1);
  }
 }
 if(line == NULL || *line == EOF) {
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
 system->bondList = malloc(sizeof(int*)*nAtoms);
 system->atomNames = malloc(sizeof(char*)*nAtoms);
 // 1d arrays
 system->atomTypes = malloc(sizeof(REAL)*nAtoms*3);
 system->X = malloc(sizeof(REAL)*nAtoms*3);
 system->F = malloc(sizeof(REAL)*nAtoms*3);
 system->lambdas = malloc(sizeof(REAL)*nAtoms*3);
 system->protons = malloc(sizeof(REAL)*nAtoms*3);
 system->valence = malloc(sizeof(REAL)*nAtoms*3);
 // Spatial
 system->boxDim = malloc(sizeof(REAL*)*3);
 for(int i = 0; i < 3; i++) {
  system->boxDim = malloc(sizeof(REAL)*3);
 }
 // Read atom lines
 Vector* bonded = vectorCreate(sizeof(int), 10, NULL, INT);
 for(int i = 0; i < nAtoms; i++) {
  fgets(line, lineSize, f);
  int atomIndex = atoi(strtok(line, " "))-1;
  system->atomNames[atomIndex] = strdup(strtok(NULL, " "));
  system->X[atomIndex] = atof(strtok(NULL, " "));
  system->X[atomIndex+1] = atof(strtok(NULL, " "));
  system->X[atomIndex+2] = atof(strtok(NULL, " "));
  system->atomTypes[atomIndex] = atoi(strtok(NULL, " "));
  char* str = strtok(NULL, " ");
  while(str != NULL) {
   int bondedID = atoi(str);
   vectorAppend(bonded, &bondedID);
   str = strtok(NULL, " ");
  }
  system->bondList[atomIndex] = vectorCopy(bonded);
  vectorClear(bonded);
  assert(atomIndex == i);
  free(line);
  line = malloc(sizeof(char)*lineSize);
 }
 printXYZ(system);
 free(line);
};

void printXYZ(System* system) {
 assert(system != NULL);
 for(int i = 0; i < system->nAtoms; i++) {
  double x = system->X[i];
  double y = system->X[i+1];
  double z = system->X[i+2];
  printf("Atom %d Name %s Type %d R=(%lf,%lf,%lf) Bonded=[", i+1, system->atomNames[i], system->atomTypes[i], x, y, z);
  Vector* bonded = system->bondList[i];
  for(int j = 0; j < bonded->size; j++) {
   int bondedAtomID = ((int*)bonded->array)[j];
   printf("%d,", bondedAtomID);
  }
  printf("]\n\n");
 }
}

void readARC(System* system, char* structureName);
void writeXYZFile(System* system, char* outputFileName);
