// Author(s): Matthew Speranza
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../include/vector.h"

int BUFFER_FACTOR = 2; // Minimum size that works for all cases

Vector *vectorCreate(int bytesPerElement, int initialCapacity, CallbackFree cbFree, enum DataType dt) {
  Vector* vec = malloc(sizeof(Vector));
  vec->size = 0;
  vec->bufSize = initialCapacity == 0 ? 1 : initialCapacity; // zero breaks resizing
  vec->bytesPerElement = bytesPerElement;
  vec->dataType = dt;
  vec->array = malloc(initialCapacity * bytesPerElement);
  vec->callbackFree = cbFree;
  return vec;
}

Vector* vectorCopy(Vector* vec) {
  Vector* retVec = malloc(sizeof(Vector));
  retVec->size = vec->size;
  retVec->bufSize = vec->bufSize;
  retVec->callbackFree = vec->callbackFree;
  retVec->bytesPerElement = vec->bytesPerElement;
  retVec->array = malloc(vec->bufSize*vec->bytesPerElement);
  memcpy(retVec->array, vec->array, vec->bufSize*vec->bytesPerElement);
  return retVec;
}

/**
 * The vector pointer that has been passed is not freed.
 * @param vec vector whose array will be freed
 */
void vectorBackingFree(Vector* vec) {
  assert(vec != NULL);
  vectorTrim(vec); // Avoid dereferencing uninitialized data in 2d case
  if(vec->callbackFree != NULL) {
    vec->callbackFree(vec->array, vec->bufSize);
  } else {
    free(vec->array);
  }
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
  vec->size = 0;
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
  if (vec->size >= vec->bufSize) { // Indexing out of bounds after this
    vectorResize(vec);
  }
  // Support 1D and 2D arrays of these basic datatypes
  switch (vec->dataType) {
    case INT: ((int *) vec->array)[vec->size++] = *(int *) elem;
      break;
    case INT_PTR: ((int **) vec->array)[vec->size++] = *(int **) elem;
      break;
    case LONG: ((long *) vec->array)[vec->size++] = *(long *) elem;
      break;
    case LONG_PTR: ((long **) vec->array)[vec->size++] = *(long **) elem;
      break;
    case FLOAT: ((float *) vec->array)[vec->size++] = *(float *) elem;
      break;
    case FLOAT_PTR: ((float **) vec->array)[vec->size++] = *(float **) elem;
      break;
    case DOUBLE: ((double *) vec->array)[vec->size++] = *(double *) elem;
      break;
    case DOUBLE_PTR: ((double **) vec->array)[vec->size++] = *(double **) elem;
      break;
    case BOOL: ((bool *) vec->array)[vec->size++] = *(bool *) elem;
      break;
    case BOOL_PTR: ((bool **) vec->array)[vec->size++] = *(bool **) elem;
      break;
    case CHAR: ((char*) vec->array)[vec->size++] = *(char *) elem;
      break;
    case CHAR_PTR: ((char**) vec->array)[vec->size++] = *(char **) elem;
      break;
    case OTHER: ((void**) vec->array)[vec->size++] = elem;
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

/**
 * Sets all data to 0 and size to 1 (but doesn't change buffer).
 * @param vec vector to set to zeros
 */
void vectorClear(Vector* vec) {
  assert(vec != NULL);
  memset(vec->array, 0, vec->bufSize*vec->bytesPerElement);
  if(vec->array == NULL) {
    printf("memset() failed in vector.c!");
    exit(1);
  }
  vec->size = 0;
}

void vectorTrim(Vector *vec) {
  assert(vec != NULL);
  vec->bufSize = vec->size;
  if(vec->bufSize == 0) {
    vec->bufSize++;
  }
  void *arr = malloc(vec->bufSize * vec->bytesPerElement);
  memcpy(arr, vec->array, vec->bufSize * vec->bytesPerElement);;
  if (vec->callbackFree) {
    vec->callbackFree(vec->array, vec->size);
  } else {
    free(vec->array);
  }
  vec->array = arr;
}

/**
 * Takes in size so that user can choose to print size or bufSize.
 * @param vec vector to print
 * @param size which length to print
 */
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
      vectorPrint(vec, vec->size);
      vectorPrint(vec, vec->bufSize);
      printf("\n");
    }
  }
  for(int i = 0; i < 10; i++) {
    num -= i;
    vectorAppend(vec, &num);
    if(verbose) {
      vectorPrint(vec, vec->size);
      vectorPrint(vec, vec->bufSize);
      printf("\n");
    }
  }
  assert(vec->size == 12);
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
  vectorBackingFree(vec);
  printf("All tests of vector.c passed!\n");
}
