//
// Created by matthew-speranza on 8/18/24.
//

#include "../include/forceFieldReader.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../system/system.h"

enum ForceFieldParams stringToFFTermEnum(char* ffTerm) {
  int index = -1;
  for(int i = 0; i < 24; i++) { // length needs to be changed if terms change
    if(strcasecmp(ForceFieldParamsStr[i], ffTerm) == 0) {
      index = i;
      break;
    }
  }
  enum ForceFieldParams param = index;
  return param;
}

void readFFLine(Vector* vec, ForceField* ff, char* line, FILE* file) {
  assert(vec != NULL);
  char** words = vec->array;
  int size = vec->size;
  assert(size > 0);
  char* command = words[0];
  enum ForceFieldParams param = stringToFFTermEnum(command);
  switch (param) {
    case ATOM: break;
    case ANGLE: break;
    case ANGLEP: break;
    case ANGTORS: break;
    case BIOTYPE: break;
    case BOND: break;
    case CHARGE: break;
    case MULTIPOLE: break;
    case OPBEND: break;
    case STRBND: break;
    case PITORS: break;
    case IMPTORS: break;
    case STRTORS: break;
    case TORSION: break;
    case IMPROPER: break;
    case TORTORS: break;
    case UREYBRAD: break;
    case VDW: break;
    case VDW14: break;
    case VDWPR: break;
    case VDWPAIR: break;
    case POLARIZE: break;
    case RELATIVESOLV: break;
    case SOLUTE: break;
    default:
      if(strcasecmp(command, "forcefield") == 0) {

      }
      break;
  }
}

void initForceField(ForceField* ff) {
  assert(ff != NULL);
  ff->atom = vectorCreate(sizeof(Atom), 20, NULL, OTHER);
  ff->angle = vectorCreate(sizeof(Angle), 20, NULL, OTHER);
  ff->angTors = vectorCreate(sizeof(AngTors), 20, NULL, OTHER);
  ff->bioType = vectorCreate(sizeof(BioType), 20, NULL, OTHER);
  ff->bond = vectorCreate(sizeof(Bond), 20, NULL, OTHER);
  ff->multipole = vectorCreate(sizeof(Multipole), 20, NULL, OTHER);
  ff->opBend = vectorCreate(sizeof(OPBend), 20, NULL, OTHER);
  ff->strBend = vectorCreate(sizeof(StrBend), 20, NULL, OTHER);
  ff->piTors = vectorCreate(sizeof(PiTors), 20, NULL, OTHER);
  ff->impTors = vectorCreate(sizeof(ImpTors), 20, NULL, OTHER);
  ff->strTors = vectorCreate(sizeof(StrTors), 20, NULL, OTHER);
  ff->torsion = vectorCreate(sizeof(Torsion), 20, NULL, OTHER);
  ff->torsTors = vectorCreate(sizeof(TorsTors), 20, NULL, OTHER);
  ff->uRayBrad = vectorCreate(sizeof(UReyBrad), 20, NULL, OTHER);
  ff->vdw = vectorCreate(sizeof(VdW), 20, NULL, OTHER);
  ff->vdwPair = vectorCreate(sizeof(VdWPair), 20, NULL, OTHER);
  ff->polarize = vectorCreate(sizeof(Polarize), 20, NULL, OTHER);
  ff->relativeSolv = vectorCreate(sizeof(RelativeSolv), 20, NULL, OTHER);
  ff->solute = vectorCreate(sizeof(Solute), 20, NULL, OTHER);
}

void forceFieldFree(ForceField* ff) {
  vectorFree(ff->atom);
  vectorFree(ff->angle);
  vectorFree(ff->angTors);
  vectorFree(ff->bioType);
  vectorFree(ff->bond);
  vectorFree(ff->multipole);
  vectorFree(ff->opBend);
  vectorFree(ff->strBend);
  vectorFree(ff->piTors);
  vectorFree(ff->impTors);
  vectorFree(ff->strTors);
  vectorFree(ff->torsion);
  vectorFree(ff->torsTors);
  vectorFree(ff->uRayBrad);
  vectorFree(ff->vdw);
  vectorFree(ff->vdwPair);
  vectorFree(ff->polarize);
  vectorFree(ff->relativeSolv);
  vectorFree(ff->solute);
  // Free the rest
  free(ff);
}

void readForceFieldFile(ForceField* forcefield, char* forceFieldFile) {
  assert(forceFieldFile != NULL);
  int len = strlen(forceFieldFile);
  if(forceFieldFile[len-1] == '\n') {
    forceFieldFile[len-1] = 0;
  }
  FILE* file = fopen(forceFieldFile, "r");
  if(file == NULL) {
    printf("Failed to read file: %s", forceFieldFile);
    exit(1);
  }
  initForceField(forcefield);
  int lineSize = 1e3;
  char line[lineSize];
  fgets(line, lineSize, file);
  while(fgets(line, lineSize, file) != NULL) {
    if(*line == EOF) {
      break;
    }
    // Read past whitespace
    while(strcmp(&line[0], "\n") == 0 && fgets(line, lineSize, file) != NULL) {
     if(*line == EOF) {
       break;
     }
    }
    // Parse line
    char* token = strtok(line, " ");
    Vector* args = vectorCreate(sizeof(char*), 10, NULL, CHAR_PTR);
    while(token != NULL && strcasecmp("#", token) != 0 && strcasecmp("//", token) != 0) {
      vectorAppend(args, &token);
      token = strtok(NULL, " ");
    }
    readFFLine(args, forcefield, line, file);
  }
}
