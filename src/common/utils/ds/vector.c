#include "../../include/vector.h"
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int BUFFER_FACTOR = 2; // Minimum size that works for all cases

Vector *vectorCreate(int bytesPerElement, int initialCapacity, CallbackFree cbFree, enum DataType dt) {
  Vector* vec = malloc(sizeof(Vector));
  vec->nextIndex = 0;
  vec->bufSize = initialCapacity == 0 ? 1 : initialCapacity; // zero breaks resizing
  vec->bytesPerElement = bytesPerElement;
  vec->dataType = dt;
  vec->array = malloc(initialCapacity * bytesPerElement);
  vec->callbackFree = cbFree;
  return vec;
}

void vectorFree(Vector* vec) {
  assert(vec != NULL);
  if(vec->callbackFree != NULL) {
    vec->callbackFree(vec->array);
  } else {
    free(vec->array);
  }
  free(vec);
}

/**
 * @param bytesPerElement Number of bytes per element of underlying data
 * @param newCapacity Capacity of backing array
 * @param cbFree Method called to free memory
 * @param initArray Inital array that gets reallocated. This pointer is freed.
 * @return pointer to vector
 */
Vector *vectorFromArray(int bytesPerElement, int newCapacity, CallbackFree cbFree, void *initArray) {
  Vector *vec = malloc(sizeof(Vector));
  vec->nextIndex = 0;
  vec->bufSize = newCapacity;
  vec->bytesPerElement = bytesPerElement;
  vec->array = reallocarray(initArray, newCapacity, bytesPerElement);
  if(vec->array != NULL) {
    free(initArray);
  } else {
    printf("Reallocation of array from init array failed in vector.c.");
    exit(1);
  }
  vec->callbackFree = cbFree;
  return vec;
}

void vectorAppend(Vector *vec, void *elem) {
  assert(vec != NULL);
  if (vec->nextIndex >= vec->bufSize) { // Indexing out of bounds after this
    vectorResize(vec);
  }
  switch (vec->dataType) {
    case INT: ((int *) vec->array)[vec->nextIndex++] = *(int *) elem;
      break;
    case LONG: ((long *) vec->array)[vec->nextIndex++] = *(long *) elem;
      break;
    case FLOAT: ((float *) vec->array)[vec->nextIndex++] = *(float *) elem;
      break;
    case DOUBLE: ((double *) vec->array)[vec->nextIndex++] = *(double *) elem;
      break;
    case BOOL: ((bool *) vec->array)[vec->nextIndex++] = *(bool *) elem;
      break;
    default:
      // This should never happen, but its an error either way.
      printf("DataType match failed in vector.c");
      exit(1);
  }
}

void vectorResize(Vector *vec) {
  assert(vec != NULL);
  vec->bufSize = BUFFER_FACTOR * vec->bufSize;
  void *arr = reallocarray(vec->array, vec->bufSize, vec->bytesPerElement);
  if (arr != NULL) {
    vec->array = arr;
  } else {
    free(vec->array);
    printf("Failed realloc in vector.c");
    exit(1);
  }
}

void vectorTrim(Vector *vec) {
  assert(vec != NULL);
  vec->bufSize = vec->nextIndex;
  void *arr = malloc(vec->bufSize * vec->bytesPerElement);
  memcpy(arr, vec->array, vec->bufSize * vec->bytesPerElement);;
  if (vec->callbackFree) {
    vec->callbackFree(vec->array);
  } else {
    free(vec->array);
  }
  vec->array = arr;
}

void vectorPrint(Vector *vec, int size) {
  assert(vec != NULL);
  printf("[");
  for (int i = 0; i < size; i++) {
    switch (vec->dataType) {
      case INT: printf("%5d,", ((int *) vec->array)[i]);
        break;
      case LONG: printf("%5ld,", ((long *) vec->array)[i]);
        break;
      case FLOAT: printf("%7.2f,", ((float *) vec->array)[i]);
        break;
      case DOUBLE: printf("%7.2f,", ((double *) vec->array)[i]);
        break;
      case BOOL: printf("%1d,", ((bool *) vec->array)[i]);
        break;
      default:
        // This should never happen, but its an error either way.
        printf("Failed to match for vector printing in vector.c");
        exit(1);
    }
  }
  printf("]\n");
}

////////////////////////////////////////////// TESTS

void vectorTest(bool verbose) {
  Vector *vec = vectorCreate(sizeof(int), 1, NULL, INT);
  int num = 10;
  for(int i = 0; i < 2; i++) {
    vectorAppend(vec, &num);
    if(verbose) {
      vectorPrint(vec, vec->nextIndex);
      vectorPrint(vec, vec->bufSize);
      printf("\n");
    }
  }
  for(int i = 0; i < 10; i++) {
    num -= i;
    vectorAppend(vec, &num);
    if(verbose) {
      vectorPrint(vec, vec->nextIndex);
      vectorPrint(vec, vec->bufSize);
      printf("\n");
    }
  }
  assert(vec->nextIndex == 12);
  if(verbose) {
    printf("Full vector size: %3d\n", vec->bufSize);
    printf("Full vector (unset mem is at end): ");
    vectorPrint(vec, vec->bufSize);
  }
  vectorTrim(vec);
  assert(vec->bufSize == 12);
  if(verbose){
    printf("Trimmed vector: ");
    vectorPrint(vec, vec->bufSize);
  }
  vectorFree(vec);
  printf("All tests of vector.c passed!\n");
}
