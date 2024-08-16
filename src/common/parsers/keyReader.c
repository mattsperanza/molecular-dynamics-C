// Author(s): Matthew Speranza
#include <stdio.h>
#include <stdlib.h>
#include "../include/keyReader.h"

#include <assert.h>
#include <string.h>
char* MD_C_Keywords[24] =
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
 "forcefield",
 "parameters",
 "params",
 "patch"
};

void handleArgs(Vector* args) {
 assert(args != NULL);
 char** words = args->array;
 char* command = words[0];
}

void readKeyFile(System* system, char* keyFile) {
 FILE* file = fopen(keyFile, "r");
 system->keyFileName = keyFile;
 int lineSize = 1e3;
 char* line = malloc(sizeof(char)*lineSize);
 fgets(line, lineSize, file);
 while(line != NULL) {
  // Read past whitespace
  while(strcmp(&line[0], "\n") == 0) {
   fgets(line, lineSize, file);
   if(line == NULL || *line == EOF) {
    printf("Failed to read from file: %s", keyFile);
    exit(1);
   }
  }
  // Parse line
  char* str = strtok(line, " ");
  Vector* args = vectorCreate(sizeof(char*), 3, NULL, CHAR_PTR);
  while(str != NULL) {
   vectorAppend(args, &str);
   str = strtok(NULL, " ");
  }
  handleArgs(args);
  // Read new line
  vectorFree(args); // No deep free needed
  free(line); // This frees everything
  line = malloc(sizeof(char)*lineSize);
  line = fgets(line, lineSize, file);
 }
 free(line);
};
