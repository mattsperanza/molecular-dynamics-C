// Author(s): Matthew Speranza
#include <stdio.h>
#include <stdlib.h>
#include "../include/keyReader.h"

#include <assert.h>
#include <string.h>
char* MD_C_Keywords[25] =
 {"verbose",
 "dt", "dtNano", "dtAtto",
 "steps",
 "printThermoEvery",
 "printRestartEvery",
 "temp", "temperature",
 "cutoff",
 "buffer",
 "a-axis","b-axis","c-axis","boxDim",
 "pmeAlpha",
 "pmeBeta",
 "pmeOrder",
 "pmeGridCount",
 "polerization",
 "forcefield", "parameters", "params",
 "patch",
 "printArchiveEvery"
};

void handleArgs(Vector* args, System* system) {
 assert(args != NULL);
 char** words = args->array;
 int size = args->size;
 char* command = strtok(words[0], "\n"); // rm newline
 if (strcasecmp(MD_C_Keywords[0], command) == 0) {
  // verbose setting
  system->verbose = true;
 } else if (strcasecmp(MD_C_Keywords[1], command) == 0 || strcasecmp(MD_C_Keywords[2], command) == 0) {
  // timestep in nanoseconds
  if(size != 2) {
   printf("Incorrect args for dt or dtFemto!");
   exit(1);
  }
  system->dtFemto = atof(words[1]);
  system->dtAtto = (int) system->dtFemto * 1e3; // truncate
 } else if (strcasecmp(MD_C_Keywords[3], command) == 0) {
  // timestep in attoseconds
  if(size != 2) {
   printf("Incorrect args for dtAtto!");
   exit(1);
  }
  system->dtAtto = atof(words[1]);
  system->dtFemto = ((REAL) system->dtAtto) / 1e3;
 } else if (strcasecmp(MD_C_Keywords[4], command) == 0) {
  // steps
  if(size != 2) {
   printf("Incorrect args for steps!");
   exit(1);
  }
  system->steps = atol(words[1]);
 } else if (strcasecmp(MD_C_Keywords[5], command) == 0) {
  // printThermoEvery
  if(size != 2) {
   printf("Incorrect args for printThermoEvery!");
   exit(1);
  }
  system->printThermoEvery = atol(words[1]);
 } else if (strcasecmp(MD_C_Keywords[6], command) == 0) {
  // printRestartEvery
  if(size != 2) {
   printf("Incorrect args for printRestartEvery!");
   exit(1);
  }
  system->printRestartEvery = atol(words[1]);
 } else if (strcasecmp(MD_C_Keywords[7], command) == 0 || strcasecmp(MD_C_Keywords[8], command) == 0) {
  // temp || temperature
  if(size != 2) {
   printf("Incorrect args for temp or temperature!");
   exit(1);
  }
  system->temperature = atof(words[1]);
 } else if (strcasecmp(MD_C_Keywords[9], command) == 0) {
  // cutoff
  if(size != 2) {
   printf("Incorrect args for cutoff!");
   exit(1);
  }
  system->realspaceCutoff = atof(words[1]);
 } else if (strcasecmp(MD_C_Keywords[10], command) == 0) {
  // buffer
  if(size != 2) {
   printf("Incorrect args for buffer!");
   exit(1);
  }
  system->realspaceBuffer = atof(words[1]);
 } else if (strcasecmp(MD_C_Keywords[11], command) == 0 || strcasecmp(MD_C_Keywords[12], command) == 0
  || strcasecmp(MD_C_Keywords[13], command) == 0) {
  // a-axis || b-axis || c-axis -- expect 3x3 REAL matrix init in file parser
  if(size == 2) { // one dim given sets all axis
   for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
     system->boxDim[i][j] = 0.0f;
    }
   }
   system->boxDim[0][0] = atof(words[1]);
   system->boxDim[1][0] = atof(words[1]);
   system->boxDim[2][0] = atof(words[1]);
  } else if(size == 4) { // three dim given
   if(strcasecmp(MD_C_Keywords[11], command) == 0) { // a-axis
    system->boxDim[0][0] = atof(words[1]); // x
    system->boxDim[0][1] = atof(words[2]); // y
    system->boxDim[0][2] = atof(words[3]); // z
   } else if(strcasecmp(MD_C_Keywords[12], command) == 0) { // b-axis
    system->boxDim[1][0] = atof(words[1]);
    system->boxDim[1][1] = atof(words[2]);
    system->boxDim[1][2] = atof(words[3]);
   } else if(strcasecmp(MD_C_Keywords[13], command) == 0) { // c-axis
    system->boxDim[2][0] = atof(words[1]);
    system->boxDim[2][1] = atof(words[2]);
    system->boxDim[2][2] = atof(words[3]);
   }
  } else {
   printf("Incorrect args for a-axis, b-axis, and c-axis!");
   exit(1);
  }
 } else if (strcasecmp(MD_C_Keywords[14], command) == 0) {
  // boxDim
  if(size != 4) {
   printf("Incorrect args for boxDim!");
   exit(1);
  }
  for(int i = 0; i < 3; i++) {
   for(int j = 0; j < 3; j++) {
    system->boxDim[i][j] = 0.0f;
   }
   system->boxDim[0][0] = atof(words[1]);
   system->boxDim[1][1] = atof(words[2]);
   system->boxDim[2][2] = atof(words[3]);
  }
 } else if (strcasecmp(MD_C_Keywords[15], command) == 0) {
  // pme alpha
  if(size != 2) {
   printf("Incorrect args for ewaldAlpha!");
   exit(1);
  }
  system->ewaldAlpha = atof(words[1]);
 } else if (strcasecmp(MD_C_Keywords[16], command) == 0) {
  // pme beta
  if(size != 2) {
   printf("Incorrect args for ewaldBeta!");
   exit(1);
  }
  system->ewaldBeta = atof(words[1]);
 } else if (strcasecmp(MD_C_Keywords[17], command) == 0) {
  // pme order
  if(size != 2) {
   printf("Incorrect args for pmeOrder!");
   exit(1);
  }
  system->ewaldOrder = atoi(words[1]);
 } else if (strcasecmp(MD_C_Keywords[18], command) == 0) {
  // pme grid count -- expected to be malloc'ed in struct file reader
  if(size == 2) {
   system->pmeGridspace[0] = atoi(words[1]);
   system->pmeGridspace[1] = atoi(words[1]);
   system->pmeGridspace[2] = atoi(words[1]);
  } else if (size == 4) {
   system->pmeGridspace[0] = atoi(words[1]);
   system->pmeGridspace[1] = atoi(words[2]);
   system->pmeGridspace[2] = atoi(words[3]);
  } else {
   printf("Incorrect args for pmeGridCount!");
   exit(1);
  }
 } else if (strcasecmp(MD_C_Keywords[19], command) == 0) {
  // polarization
  if(size != 2) {
   printf("Incorrect args for polarization!");
   exit(1);
  }
  if(strcasecmp("NONE", words[1]) == 0) {
   system->polarization = NONE;
  } else if (strcasecmp("DIRECT", words[1]) == 0) {
   system->polarization = DIRECT;
  } else if (strcasecmp("MUTUAL", words[1]) == 0) {
   system->polarization = MUTUAL;
  }
 } else if (strcasecmp(MD_C_Keywords[20], command) == 0 || strcasecmp(MD_C_Keywords[21], command) == 0
  || strcasecmp(MD_C_Keywords[22], command) == 0) {
  // forcefield || parameters || params
  if(size != 2) {
   printf("Incorrect args for forcefield/parameters!");
   exit(1);
  }
  system->forceFieldFile = strdup(words[1]);
 } else if (strcasecmp(MD_C_Keywords[23], command) == 0) {
  // patch -- vector created in struct file reader
  if(size != 2) {
   printf("Incorrect args for patch!");
   exit(1);
  }
  vectorAppend(system->patchFiles, strdup(words[1]));
 } else if (strcasecmp(MD_C_Keywords[24], command) == 0) {
  // printArchiveEvery
  if(size != 2) {
   printf("Incorrect args for printArchiveEvery!");
   exit(1);
  }
  system->printArchiveEvery = atol(words[1]);
 }
}

void readKeyFile(System* system, char* keyFile) {
 FILE* file = fopen(keyFile, "r");
 system->keyFileName = keyFile;
 int lineSize = 1e3;
 char line[lineSize];
 fgets(line, lineSize, file);
 while(fgets(line, lineSize, file) != NULL) {
  // Read past whitespace
  while(strcmp(&line[0], "\n") == 0 && fgets(line, lineSize, file) != NULL) {
   if(*line == EOF) {
    printf("Failed to read from file: %s", keyFile);
    exit(1);
   }
  }
  // Parse line
  char* str = strtok(line, " ");
  Vector* args = vectorCreate(sizeof(char*), 3, NULL, CHAR_PTR);
  // Ignore comments
  while(str != NULL && strcasecmp("#", str) != 0 && strcasecmp("//", str) != 0) {
   vectorAppend(args, &str);
   str = strtok(NULL, " ");
  }
  handleArgs(args, system);
  // Read new line
  vectorFree(args); // No deep free needed
  fgets(line, lineSize, file);
 }
 fclose(file);
};
