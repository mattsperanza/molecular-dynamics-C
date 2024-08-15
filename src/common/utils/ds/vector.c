#include "../../include/vector.h"
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float BUFFER_FACTOR = 1.5f;

Vector *vectorCreate(int bytesPerElement, int initialCapacity, CallbackFree cbFree, enum DataType dt) {
  Vector *vec = malloc(sizeof(Vector));
  vec->nextIndex = 0;
  vec->bufSize = initialCapacity;
  vec->bytesPerElement = bytesPerElement;
  vec->dataType = dt;
  vec->array = malloc(initialCapacity * bytesPerElement);
  vec->callbackFree = cbFree;
  return vec;
}

void vectorFree(Vector *vec) {
  assert(vec != NULL);
  if (vec->callbackFree) {
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
  vec->array = realloc(initArray, newCapacity * bytesPerElement);
  if(vec->array) {
    free(initArray);
  }
  vec->callbackFree = cbFree;
  return vec;
}

void vectorAppend(Vector *vec, void *elem) {
  assert(vec != NULL);
  if (vec->nextIndex + 1 >= vec->bufSize) {
    // Indexing will go out of bounds after this
    vec->bufSize = (int) (BUFFER_FACTOR * vec->bufSize);
    void* arr = realloc(vec->array, vec->bufSize * vec->bytesPerElement);
    if(arr != NULL) {
      vec->array = arr;
    } else {
      free(vec->array);
      printf("Failed to realloc in vector.c!");
      exit(1);
    }
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
  vec->bufSize = (int) (BUFFER_FACTOR * vec->bufSize * vec->bytesPerElement);
  void *arr = realloc(vec->array, vec->bufSize);
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
  vec->bufSize = vec->nextIndex + 1;
  void *arr = malloc(vec->bufSize * vec->bytesPerElement);
  memcpy(arr, vec->array, vec->bufSize * vec->bytesPerElement);;
  if (vec->callbackFree) {
    vec->callbackFree(vec->array);
  } else {
    free(vec->array);
  }
  vec->array = arr;
}

void vectorPrint(Vector *vec) {
  assert(vec != NULL);
  printf("[");
  for (int i = 0; i < vec->nextIndex; i++) {
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

void vectorTest() {
  Vector *vec = vectorCreate(sizeof(int), 10, NULL, INT);
  int num = 10;
  vectorAppend(vec, &num);
  vectorPrint(vec);
  vectorAppend(vec, &num);
  for(int i = 0; i < 15; i++) {
    num -= i;
    vectorAppend(vec, &num);
    vectorPrint(vec);
  }
  vectorFree(vec);
  printf("All tests of vector.c passed!\n");
}
