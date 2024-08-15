#include <assert.h>
#include <stddef.h>
#include <stdio.h>

#include "../include/vector.h"

void printIntegerVector(Vector* vec);

void vectorTest() {
  Vector* vec = vectorCreate(sizeof(int), 10, NULL, INT);
  int num = 10;
  vectorAppend(vec, &num);
  printIntegerVector(vec);
}

void printIntegerVector(Vector* vec) {
  assert(vec != NULL);
  printf("[");
  for(int i = 0; i < vec->index+1; i++) {
    printf("%3d,", ((int*) vec->array)[i]);
  }
  printf("]\n");
}