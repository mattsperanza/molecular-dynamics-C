// Author(s): Matthew Speranza
#include "../include/forceFieldReader.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

int compareInt(const void* a, const void* b) {
  return *(int*)a - *(int*)b;
}

Atom* atomLine(char** words, int size) {
  Atom* atom = malloc(sizeof(Atom));
  if(atom == NULL) {
    printf("Couldn't read atom line.");
    exit(1);
  }
  if(size < 8) {
    printf("Incorrect number of arguments for atom line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  atom->type = atoi(words[1]);
  atom->aClass = atoi(words[2]);
  for(int i = 0; i < 3; i++) {
    atom->name[i] = words[3][i];
  }
  atom->name[3] = '\0';
  int end = size-3;
  int count = 0;
  for(int i = 4; i < end; i++) {
    int len = strlen(words[i]);
    for(int j = 0; j < len; j++) {
      atom->environment[count++] = words[i][j];
    }
  }
  atom->environment[count] = '\0';
  atom->atomicNum = atoi(words[end++]);
  atom->atomicMass = atof(words[end++]);
  atom->valence = atoi(words[end]);
  return atom;
}

Angle* angleLine(char** words, int size) {
  Angle* angle = malloc(sizeof(Angle));
  if(angle == NULL) {
    printf("Couldn't allocate angle line.");
    exit(1);
  }
  if (size < 5) {
    printf("Couldn't read angle line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  for(int i = 0; i < 3; i++) {
    angle->aClasses[i] = atoi(words[i+1]);
  }
  qsort(angle->aClasses, 3, sizeof(int), compareInt);
  angle->forceConstant = atof(words[4]);
  assert(size-5 > 0 && size-5 < 4);
  for(int i = 0; i < size-5; i++) {
    angle->angle[i] = atof(words[i+5]);
  }
  return angle;
}

AngTors* angtorsLine(char** words, int size) {
  AngTors* angtors = malloc(sizeof(Angle));
  if(angtors == NULL) {
    printf("Couldn't allocate angle line.");
    exit(1);
  }
  if (size != 11) {
    printf("Couldn't read angle line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  for(int i = 0; i < 4; i++) {
    angtors->aClasses[i] = atoi(words[i+1]);
  }
  qsort(angtors->aClasses, 4, sizeof(int), compareInt);
  for(int i = 0; i < 6; i++) {
    angtors->forceConstants[i] = atof(words[i+5]);
  }
  return angtors;
}

BioType* biotypeLine(char** words, int size) {
  BioType* biotype = malloc(sizeof(BioType));
  if (biotype == NULL) {
    printf("Couldn't allocate biotype line.");
    exit(1);
  }
  if(size < 5) {
    printf("Too few arguments for biotype line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  biotype->index = atoi(words[1]);
  int len = strlen(words[2]);
  assert(len <= 4);
  for(int i = 0; i < len; i++) {
    biotype->atomName[i] = words[2][i];
  }
  biotype->atomName[len] = '\0';
  int end = size-1; // len-(#args after moleculeName)
  int count = 0;
  for(int i = 3; i < end; i++) {
    int len = strlen(words[i]);
    for(int j = 0; j < len; j++) {
      biotype->moleculeName[count++] = words[i][j];
    }
  }
  biotype->moleculeName[count] = '\0';
  biotype->atomType = atof(words[end]);
  return biotype;
}

Bond* bondLine(char** words, int size) {
  Bond* bond = malloc(sizeof(Bond));
  if(bond == NULL) {
    printf("Couldn't allocate bond line.");
    exit(1);
  }
  if(size != 5) {
    printf("Too few arguments for bond line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  bond->atomClasses[0] = atoi(words[1]);
  bond->atomClasses[1] = atoi(words[2]);
  qsort(bond->atomClasses, 2, sizeof(int), compareInt);
  bond->forceConstant = atof(words[3]);
  bond->distance = atof(words[4]);
  return bond;
}

Multipole* multipoleLines(char** words, int size, char* line, FILE* file) {
  Multipole* mpole = malloc(sizeof(Multipole));
  if (mpole == NULL) {
    printf("Couldn't allocate multipole line.");
    exit(1);
  }
  if(strcasecmp(words[0], "charge") == 0) {
    if(size != 3) {
      printf("Couldn't parse charge line: ");
      for(int i = 0; i < size; i++) {
        printf("%s ", words[i]);
      }
      printf("\n");
      exit(1);
    }
    mpole->frameAtomTypes[0] = atoi(words[1]);
    mpole->multipole[0] = atof(words[2]);
  } else {
    if(size < 5) {
      printf("Couldn't parse multipole line: ");
      for(int i = 0; i < size; i++) {
        printf("%s ", words[i]);
      }
      printf("\n");
      exit(1);
    }
    mpole->frameAtomTypes[0] = atoi(words[1]);
    mpole->frameAtomTypes[1] = atoi(words[2]);
    mpole->frameAtomTypes[2] = atoi(words[3]);
    int next = 4;
    if(size != 5) {
      mpole->frameAtomTypes[3] = atoi(words[next++]);
    }
    if(next == 4) { // 3 atom classes
      if(mpole->frameAtomTypes[1] < 0 || mpole->frameAtomTypes[2] < 0) {
        mpole->frameDef = BISECTOR;
      } else {
        mpole->frameDef = ZTHENX;
      }
    } else if(next == 5) { // 4 atom classes
      if(mpole->frameAtomTypes[2] < 0 || mpole->frameAtomTypes[3] < 0) {
        mpole->frameDef = ZTHENBISECTOR;
      } else {
        mpole->frameDef = THREEFOLD;
      }
    }
    if(mpole->frameDef == NONE) {
      printf("Frame definition not defined for multipole line: ");
      for(int i = 0; i < size; i++) {
        printf("%s ", words[i]);
      }
      printf("\n");
      exit(1);
    }
    // mpole = [q, dx, dy, dz, qxx, qyy, qzz, 2*qxy, 2*qxz, 2*qyz]
    mpole->multipole[0] = atof(words[next++]); // charge
    fgets(line, 1e3, file);
    mpole->multipole[1] = atof(strtok(line, " ")); // dx
    mpole->multipole[2] = atof(strtok(NULL, " ")); // dy
    mpole->multipole[3] = atof(strtok(NULL, " ")); // dz
    fgets(line, 1e3, file);
    mpole->multipole[4] = atof(strtok(line, " "))/3; // qxx
    fgets(line, 1e3, file);
    mpole->multipole[7] = 2*atof(strtok(line, " "))/3; // qxy
    mpole->multipole[5] = atof(strtok(NULL, " "))/3; // qyy
    fgets(line, 1e3, file);
    mpole->multipole[8] = 2*atof(strtok(line, " "))/3; // 2*qxz
    mpole->multipole[9] = 2*atof(strtok(NULL, " "))/3; // 2*qyz
    mpole->multipole[6] = atof(strtok(NULL, " "))/3; // qzz
  }
  return mpole;
}

OPBend* opbendLine(char** words, int size) {
  OPBend* opbend = malloc(sizeof(OPBend));
  if(opbend == NULL) {
    printf("Couldn't allocate opbend line.");
    exit(1);
  }
  if(size != 6) {
    printf("Couldn't parse opbend line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  for(int i = 0; i < 4; i++) {
    opbend->atomClasses[i] = atoi(words[i+1]);
  }
  qsort(opbend->atomClasses, 4, sizeof(int), compareInt);
  opbend->forceConstant = atof(words[5]);
  return opbend;
}

StrBend* strbendLine(char** words, int size) {
  StrBend* strbend = malloc(sizeof(StrBend));
  if(strbend == NULL) {
    printf("Couldn't allocate strbend line.");
    exit(1);
  }
  if(size != 6) {
    printf("Couldn't parse strbend line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  for(int i = 0; i < 3; i++) {
    strbend->atomClasses[i] = atoi(words[i+1]);
  }
  qsort(strbend->atomClasses, 3, sizeof(int), compareInt);
  for(int i = 0; i < 2; i++) {
    strbend->forceConstants[i] = atof(words[i+4]);
  }
  return strbend;
}

PiTors* pitorsLine(char** words, int size) {
  PiTors* pitors = malloc(sizeof(PiTors));
  if(pitors == NULL) {
    printf("Couldn't allocate pitors line.");
    exit(1);
  }
  if(size != 4) {
    printf("Couldn't parse pitors line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  for(int i = 0; i < 2; i++) {
    pitors->atomClasses[i] = atoi(words[i+1]);
  }
  qsort(pitors->atomClasses, 2, sizeof(int), compareInt);
  pitors->forceConstant = atof(words[3]);
  return pitors;
}

ImpTors* imptorsLine(char** words, int size) {
  ImpTors* imptors = malloc(sizeof(ImpTors));
  if(imptors == NULL) {
    printf("Couldn't allocate imptors line.");
    exit(1);
  }
  if(size != 7) {
    printf("Couldn't parse imptors line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  for(int i = 0; i < 4; i++) {
    imptors->atomClasses[i] = atoi(words[i+1]);
  }
  qsort(imptors->atomClasses, 4, sizeof(int), compareInt);
  imptors->forceConstant = atof(words[5]);
  imptors->phase = atof(words[6]);
  imptors->periodicity = atof(words[7]);
  return imptors;
}

StrTors* strtorsLine(char** words, int size) {
  StrTors* strtors = malloc(sizeof(StrTors));
  if(strtors == NULL) {
    printf("Couldn't allocate strtors line.");
    exit(1);
  }
  if(size != 14) {
    printf("Couldn't parse strtors line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  for(int i = 0; i < 4; i++) {
    strtors->atomClasses[i] = atoi(words[i+1]);
  }
  qsort(strtors->atomClasses, 4, sizeof(int), compareInt);
  for(int i = 0; i < 9; i++) {
    strtors->forceConstants[i] = atof(words[i+5]);
  }
  return strtors;
}

Torsion* torsionLine(char** words, int size, enum TorsionMode param) {
  Torsion* torsion = malloc(sizeof(Torsion));
  if(torsion == NULL) {
    printf("Couldn't allocate torsion line.");
    exit(1);
  }
  if(size < 5) {
    printf("Couldn't parse torsion line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  for(int i = 0; i < 4; i++) {
    torsion->atomClasses[i] = atoi(words[i+1]);
  }
  qsort(torsion->atomClasses, 4, sizeof(int), compareInt);
  torsion->terms = 0;
  for(int i = 0; i < 3; i++) {
    torsion->amplitude[i] = atof(words[3*i+5]);
    if(torsion->amplitude[i] != 0.0) {torsion->terms++;}
    torsion->phase[i] = atof(words[3*i+6]);
    torsion->periodicity[i] = atoi(words[3*i+7]);
  }
  torsion->torsionMode = param;
  return torsion;
}

TorTors* tortorsLines(char** words, int size, char* line, FILE* file) {
  TorTors* tortors = malloc(sizeof(TorTors));
  //TODO: not implemented
  return tortors;
}

UReyBrad* uraybradLine(char** words, int size) {
  UReyBrad* uraybrad = malloc(sizeof(UReyBrad));
  if(uraybrad == NULL) {
    printf("Couldn't allocate uraybrad line.");
    exit(1);
  }
  if(size != 6) {
    printf("Couldn't parse uraybrad line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  for(int i = 0; i < 3; i++) {
    uraybrad->atomClasses[i] = atoi(words[i+1]);
  }
  qsort(uraybrad->atomClasses, 3, sizeof(int), compareInt);
  uraybrad->forceConstant = atof(words[4]);
  uraybrad->distance = atof(words[5]);
  return uraybrad;
}

VdW* vdwLine(char** words, int size, enum VdWType param) {
  VdW* vdw = malloc(sizeof(VdW));
  if(vdw == NULL) {
    printf("Couldn't allocate memory for vdw line!");
    exit(1);
  }
  if(size == 4 || size == 5) {
    vdw->atomClass = atoi(words[1]);
    vdw->radius = atof(words[2]);
    vdw->wellDepth = atof(words[3]);
    vdw->vdwType = param;
    if(size == 5){ vdw->reductionFactor = atof(words[4]); }
  } else {
    printf("Can't parse vdw line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  return vdw;
}

VdWPair* vdwpairLine(char** words, int size) {
  VdWPair* vdwpair = malloc(sizeof(VdWPair));
  if(vdwpair == NULL) {
    printf("Couldn't allocate memory for vdwpair line!");
    exit(1);
  }
  if(size != 5) {
    printf("Couldn't parse vdwpair line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  for(int i = 0; i < 2; i++) {
    vdwpair->atomClasses[i] = atoi(words[i+1]);
  }
  qsort(vdwpair->atomClasses, 2, sizeof(int), compareInt);
  vdwpair->radius = atof(words[3]);
  vdwpair->wellDepth = atof(words[4]);
  return vdwpair;
}

Polarize* polarizeLine(char** words, int size) {
  Polarize* polarize = malloc(sizeof(Polarize));
  if(polarize == NULL) {
    printf("Couldn't allocate memory for polarize line!");
    exit(1);
  }
  if(size < 4) {
    printf("Failed to parse polarize line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
  }
  polarize->atomType = atoi(words[1]);
  for(int i = 0; i < 3; i++) {
    polarize->polarizabilityTensor[i][i] = atof(words[2]);
  }
  polarize->thole = atof(words[3]);
  for(int i = 4; i < size; i++) {
    polarize->polarizationGroup[i-4] = atoi(words[i]);
  }
  qsort(polarize->polarizationGroup, size-4, sizeof(int), compareInt);
  return polarize;
}

RelativeSolv* relativesolvLine(char** words, int size) {
  RelativeSolv* relativesolv = malloc(sizeof(RelativeSolv));
  if(relativesolv == NULL) {
    printf("Couldn't allocate memory for relativesolv line!");
    exit(1);
  }
  return relativesolv;
}

Solute* soluteLine(char** words, int size) {
  Solute* solute = malloc(sizeof(Solute));
  if(solute == NULL) {
    printf("Could not allocate memory for solute line!");
    exit(1);
  }
  if(size != 6) {
    printf("Couldn't parse solute line: ");
    for(int i = 0; i < size; i++) {
      printf("%s ", words[i]);
    }
    printf("\n");
    exit(1);
  }
  solute->atomType = atoi(words[1]);
  for(int i = 0; i < 3; i++) {
    solute->diameters[i] = atof(words[i+2]);
  }
  solute->sneck = atof(words[5]);
  return solute;
}

void readFFLine(Vector* vec, ForceField* ff, char* line, FILE* file) {
  assert(vec != NULL);
  assert(vec->size > 0);
  char** words = vec->array;
  char* command = words[0];
  enum ForceFieldParams param = stringToFFTermEnum(command);
  switch (param) {
    case ATOM: vectorAppend(&ff->atom, atomLine(words, vec->size));
      break;
    case ANGLE: vectorAppend(&ff->angle, angleLine(words, vec->size));
      break;
    case ANGLEP: vectorAppend(&ff->angle, angleLine(words, vec->size));
      break;
    case ANGTORS: vectorAppend(&ff->angTors, angtorsLine(words, vec->size));
      break;
    case BIOTYPE: vectorAppend(&ff->bioType, biotypeLine(words, vec->size));
      break;
    case BOND: vectorAppend(&ff->bond, bondLine(words, vec->size));
      break;
    case CHARGE:
      vectorAppend(&ff->multipole, multipoleLines(words, vec->size, line, file));
      ff->name = CHARMM;
      break;
    case MULTIPOLE:
      vectorAppend(&ff->multipole, multipoleLines(words, vec->size, line, file));
      ff->name = AMOEBA;
      break;
    case OPBEND: vectorAppend(&ff->opBend, opbendLine(words, vec->size));
      break;
    case STRBND: vectorAppend(&ff->strBend, strbendLine(words, vec->size));
      break;
    case PITORS: vectorAppend(&ff->piTors, pitorsLine(words, vec->size));
      break;
    case IMPTORS: vectorAppend(&ff->impTors, imptorsLine(words, vec->size));
      break;
    case STRTORS: vectorAppend(&ff->strTors, strtorsLine(words, vec->size));
      break;
    case TORSION: vectorAppend(&ff->torsion, torsionLine(words, vec->size, TORS_NORMAL));
      break;
    case IMPROPER: vectorAppend(&ff->torsion, torsionLine(words, vec->size, TORS_IMPROPER));
      break;
    case TORTORS: vectorAppend(&ff->torTors, tortorsLines(words, vec->size, line, file));
      break;
    case UREYBRAD: vectorAppend(&ff->uRayBrad, uraybradLine(words, vec->size));
      break;
    case VDW: vectorAppend(&ff->vdw, vdwLine(words, vec->size, VDW_NORMAL));
      break;
    case VDW14: vectorAppend(&ff->vdw, vdwLine(words, vec->size, VDW_14));
      break;
    case VDWPR: vectorAppend(&ff->vdwPair, vdwpairLine(words, vec->size));
      break;
    case VDWPAIR: vectorAppend(&ff->vdwPair, vdwpairLine(words, vec->size));
      break;
    case POLARIZE: vectorAppend(&ff->polarize, polarizeLine(words, vec->size));
      break;
    case RELATIVESOLV: vectorAppend(&ff->relativeSolv, relativesolvLine(words, vec->size));
      break;
    case SOLUTE: vectorAppend(&ff->solute, soluteLine(words, vec->size));
      break;
    default:
      break;
  }
}

void initForceField(ForceField* ff) {
  assert(ff != NULL);
  ff->atom = *vectorCreate(sizeof(Atom*), 20, NULL, OTHER);
  ff->angle = *vectorCreate(sizeof(Angle*), 20, NULL, OTHER);
  ff->angTors = *vectorCreate(sizeof(AngTors*), 20, NULL, OTHER);
  ff->bioType = *vectorCreate(sizeof(BioType*), 20, NULL, OTHER);
  ff->bond = *vectorCreate(sizeof(Bond*), 20, NULL, OTHER);
  ff->multipole = *vectorCreate(sizeof(Multipole*), 20, NULL, OTHER);
  ff->opBend = *vectorCreate(sizeof(OPBend*), 20, NULL, OTHER);
  ff->strBend = *vectorCreate(sizeof(StrBend*), 20, NULL, OTHER);
  ff->piTors = *vectorCreate(sizeof(PiTors*), 20, NULL, OTHER);
  ff->impTors = *vectorCreate(sizeof(ImpTors*), 20, NULL, OTHER);
  ff->strTors = *vectorCreate(sizeof(StrTors*), 20, NULL, OTHER);
  ff->torsion = *vectorCreate(sizeof(Torsion*), 20, NULL, OTHER);
  ff->torTors = *vectorCreate(sizeof(TorTors*), 20, NULL, OTHER);
  ff->uRayBrad = *vectorCreate(sizeof(UReyBrad*), 20, NULL, OTHER);
  ff->vdw = *vectorCreate(sizeof(VdW*), 20, NULL, OTHER);
  ff->vdwPair = *vectorCreate(sizeof(VdWPair*), 20, NULL, OTHER);
  ff->polarize = *vectorCreate(sizeof(Polarize*), 20, NULL, OTHER);
  ff->relativeSolv = *vectorCreate(sizeof(RelativeSolv*), 20, NULL, OTHER);
  ff->solute = *vectorCreate(sizeof(Solute*), 20, NULL, OTHER);
}

void forceFieldFree(ForceField* ff) {
  vectorBackingFree(&ff->atom);
  vectorBackingFree(&ff->angle);
  vectorBackingFree(&ff->angTors);
  vectorBackingFree(&ff->bioType);
  vectorBackingFree(&ff->bond);
  vectorBackingFree(&ff->multipole);
  vectorBackingFree(&ff->opBend);
  vectorBackingFree(&ff->strBend);
  vectorBackingFree(&ff->piTors);
  vectorBackingFree(&ff->impTors);
  vectorBackingFree(&ff->strTors);
  vectorBackingFree(&ff->torsion);
  vectorBackingFree(&ff->torTors);
  vectorBackingFree(&ff->uRayBrad);
  vectorBackingFree(&ff->vdw);
  vectorBackingFree(&ff->vdwPair);
  vectorBackingFree(&ff->polarize);
  vectorBackingFree(&ff->relativeSolv);
  vectorBackingFree(&ff->solute);
  for(int i = 0; i < ff->nAtoms; i++) {
    free(ff->wellDepths[i]);
    free(ff->rMin[i]);
  }
  free(ff->wellDepths);
  free(ff->rMin);
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
    char* str = strtok(line, " ");
    Vector* args = vectorCreate(sizeof(char*), 3, NULL, CHAR_PTR);
    // Ignore comments & tokenize
    while(str != NULL && strcasecmp("#", &str[0]) != 0 && strcasecmp("/", &str[0]) != 0) {
      vectorAppend(args, &str);
      str = strtok(NULL, " ");
    }
    readFFLine(args, forcefield, line, file);
  }
}

void vdwParameters(ForceField* forceField, int nAtoms, int* atomClasses, Vector* neighborList) {
  printf("Calculating VdW parameters\n");
  // Assign VdW parameters to every atom
  if(forceField->wellDepths != NULL) {
    free(forceField->wellDepths);
  }
  if(forceField->rMin != NULL) {
    free(forceField->rMin);
  }
  forceField->nAtoms = nAtoms;
  forceField->wellDepths = malloc(sizeof(REAL*) * nAtoms);
  forceField->rMin = malloc(sizeof(REAL*) * nAtoms);
  forceField->reductionFactors = calloc(sizeof(REAL), nAtoms);
  for(int i = 0; i < nAtoms; i++) {
    int* list = neighborList[i].array;
    int size = neighborList[i].size;
    forceField->wellDepths[i] = malloc(sizeof(REAL) * size);
    forceField->rMin[i] = malloc(sizeof(REAL) * size);

    int iClass = atomClasses[i];
    VdW iVdW = *((VdW**)forceField->vdw.array)[iClass-1];
    REAL ri = iVdW.radius;
    REAL ri2 = ri * ri;
    REAL epsi = iVdW.wellDepth;
    REAL redi = iVdW.reductionFactor;
    forceField->reductionFactors[i] = redi;

    // Calculate VdW parameters for each atom pair in neighbor list
    for(int jj = 0; jj < size; jj++) {
      int jClass = atomClasses[list[jj]];
      // Replace VdW param with VdWPair params
      int vdwPairSize = forceField->vdwPair.size;
      VdW jVdW = *((VdW**)forceField->vdw.array)[jClass-1];
      REAL rj = jVdW.radius;
      REAL rj2 = rj * rj;
      REAL epsj = jVdW.wellDepth;
      for(int i = 0; i < vdwPairSize; i++) {
        VdWPair* vdwPair = ((VdWPair**)forceField->vdwPair.array)[i];
        if(iClass == vdwPair->atomClasses[0] || iClass == vdwPair->atomClasses[1]) {
          if(jClass == vdwPair->atomClasses[0] || jClass == vdwPair->atomClasses[1]) {
            ri = vdwPair->radius;
            ri2 = ri * ri;
            epsj = vdwPair->wellDepth;
            rj = vdwPair->radius;
            rj2 = rj * rj;
            epsj = vdwPair->wellDepth;
          }
        }
      }
      if(forceField->name == AMOEBA) {
        // (HHG) rmin.ij = (ri^3 + rj^3) / (ri^2 + rj^2)
        forceField->rMin[i][jj] = (ri2*ri + rj2*rj)/(ri2 + rj2);
        // (HHG) eps.ij  = 4.0*(eps.i * eps.j) / ((sqrt(eps.i) + sqrt(eps.j))^2)
        forceField->wellDepths[i][jj] = 4.0 * epsi * epsj / ((sqrt(epsi) + sqrt(epsj)) * (sqrt(epsi) + sqrt(epsj)));
      } else if(forceField->name == CHARMM) {
        // (GEOMETRIC) rmin.ij = 2.0 * sqrt(ri) * sqrt(rj)
        forceField->rMin[i][jj] = 2.0 * sqrt(ri) * sqrt(rj);
        // (GEOMETRIC) eps.ij  = sqrt(eps.i)*sqrt(eps.j)
        forceField->wellDepths[i][jj] = sqrt(epsi) * sqrt(epsj);
      }
    }
  }
  if(forceField->name == AMOEBA) {
    forceField->vdwM = 7;
    forceField->vdwN = 14;
    forceField->vdwDelta = 0.07;
    forceField->vdwGamma = 0.12;
  } else if (forceField->name == CHARMM) {
    forceField->vdwM = 6;
    forceField->vdwN = 12;
    forceField->vdwDelta = 0.0;
    forceField->vdwGamma = 0.0;
  }
}

void assignMultipoles(ForceField* forceField, REAL** multipoles, REAL** frameDef, Vector* list12, Vector* list13,
  int* atomClasses, int nAtoms) {
  multipoles = malloc(sizeof(REAL*)*nAtoms);
  frameDef = malloc(sizeof(int*)*nAtoms); // int[5] of atom ids (with signs kept) and last is enum of def type

  for(int i = 0; i < nAtoms; i++) {
    int atomClass = atomClasses[i];
    int* bonded = list12[i].array;
    Vector mpoles = forceField->multipole;
    multipoles[i] = NULL;

    // 0 reference atoms
    bool breakFlag = false;
    for(int j = 0; j < mpoles.size; j++) {
      int* frameAtoms = ((Multipole**)mpoles.array)[j]->frameAtomTypes;
      if(frameAtoms[0] == atomClass && frameAtoms[1] == 0 &&
        frameAtoms[2] == 0 && frameAtoms[3] == 0) {
        multipoles[i] = ((Multipole**)mpoles.array)[j]->multipole;
        breakFlag = true;
        break;
      }
    }
    if(breakFlag) {
      continue;
    }

    // 1 reference atom
    breakFlag = false;
    for(int j = 0; j < list12[i].size; j++) {
      int atomID2 = bonded[j];
      int atomClass2 = atomClasses[atomID2];
      // Loop over multipoles with two frame atoms
      for(int k = 0; k < mpoles.size; k++) {
        int* frameAtoms = ((Multipole**)mpoles.array)[k]->frameAtomTypes;
        if(frameAtoms[0] == atomClass && frameAtoms[1] == atomClass2 &&
          frameAtoms[2] == 0 && frameAtoms[3] == 0) {
          multipoles[i] = ((Multipole**)mpoles.array)[k]->multipole;
          breakFlag = true;
          break;
        }
      }
      if(breakFlag) {
        break;
      }
    }
    if(breakFlag) {
      continue;
    }

    // 2 reference atoms
    breakFlag = false;
    for(int j = 0; j < list12[i].size; j++) {
      int atomID2 = bonded[j];
      int atomClass2 = atomClasses[atomID2];
      for(int k = 0; k < list12[i].size; k++) {
        int atomID3 = bonded[k];
        int atomClass3 = atomClasses[atomID3];
        if(atomID3 != atomID2) {
          int* frameAtoms = ((Multipole**)mpoles.array)[k]->frameAtomTypes;
          if(frameAtoms[0] == atomClass && frameAtoms[1] == atomClass2 &&
            frameAtoms[2] == atomClass3 && frameAtoms[3] == 0) {
            multipoles[i] = ((Multipole**)mpoles.array)[k]->multipole;
            breakFlag = true;
            break;
          }
        }
      }
      if(breakFlag) {
        break;
      }
    }
    if(breakFlag) {
      continue;
    }

    // 3 reference atoms
    breakFlag = false;
    for(int j = 0; j < list12[i].size; j++) {
      int atomID2 = bonded[j];
      int atomClass2 = atomClasses[atomID2];
      for(int k = 0; k < list12[i].size; k++) {
        int atomID3 = bonded[k];
        int atomClass3 = atomClasses[atomID3];
        if(atomID3 != atomID2) {
          for(int l = 0; l < list12[i].size; l++) {
            int atomID4 = bonded[l];
            int atomClass4 = atomClasses[atomID4];
            if(atomID4 != atomID3 || atomID4 != atomID2) {
              int* frameAtoms = ((Multipole**)mpoles.array)[l]->frameAtomTypes;
              if(frameAtoms[0] == atomClass && frameAtoms[1] == atomClass2 &&
                frameAtoms[2] == atomClass3 && frameAtoms[3] == atomClass4) {
                multipoles[i] = ((Multipole**)mpoles.array)[l]->multipole;
                breakFlag = true;
                break;
              }
            }
          }
        }
        if(breakFlag) {
          break;
        }
      }
      if(breakFlag) {
        break;
      }
    }
    if(breakFlag) {
      continue;
    }

    // 1-3 definition
    breakFlag = false;
    int* farSite = list13[i].array;
    for(int j = 0; j < list12[i].size; j++) {
      int atomID2 = bonded[j];
      int atomClass2 = atomClasses[atomID2];
      for(int k = 0; k < list13[i].size; k++) {
        int atomID3 = farSite[k];
        int atomClass3 = atomClasses[atomID3];
        for(int l = 0; l < mpoles.size; l++) {
          int* frameAtoms = ((Multipole**)mpoles.array)[l]->frameAtomTypes;
          if(frameAtoms[0] == atomClass && frameAtoms[1] == atomClass2 &&
            frameAtoms[2] == atomClass3 && frameAtoms[3] == 0) {
            multipoles[i] = ((Multipole**)mpoles.array)[l]->multipole;
            breakFlag = true;
            break;
          }
        }
        if(breakFlag) {
          break;
        }
      }
      if(breakFlag) {
        break;
      }
    }

    if(multipoles[i] == NULL) {
      printf("Failed to match multipole class %d for atom %d\n", atomClass, i);
      exit(1);
    }
  }

}