#include "../../include/vector.h"
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

Vector* vectorCreate(int bytesPerElement, int initialCapacity, CallbackFree cbFree, enum DataType dt) {
  Vector* vec = malloc(sizeof(Vector));
  vec->index = 0;
  vec->bufSize = initialCapacity;
  vec->bytesPerElement = bytesPerElement;
  vec->dataType = dt;
  vec->array = malloc(initialCapacity * bytesPerElement);
  vec->callbackFree = cbFree;
  return vec;
}

void freeVector(Vector* vec) {
  assert(vec != NULL);
  if(vec->callbackFree) {
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
 * @param initArray Inital array that gets reallocated
 * @return pointer to vector
 */
Vector* vectorFromArray(int bytesPerElement, int newCapacity, CallbackFree cbFree, void* initArray) {
  Vector* vec = malloc(sizeof(Vector));
  vec->index = 0;
  vec->bufSize = newCapacity;
  vec->bytesPerElement = bytesPerElement;
  vec->array = realloc(initArray, newCapacity * bytesPerElement);
  vec->callbackFree = cbFree;
  return vec;
}

void vectorAppend(Vector* vec, void* elem) {
  assert(vec != NULL);
  if(vec->index+1 >= vec->bufSize) { // Indexing will go out of bounds after this
    vec->bufSize = (int) (1.25*vec->bufSize);
    realloc(vec->array, vec->bufSize*vec->bytesPerElement);
  }
  switch (vec->dataType) {
    case INT: ((int*)vec->array)[++vec->index] = *(int*)elem;
      break;
    case LONG: ((long*)vec->array)[++vec->index] = *(long*)elem;
      break;
    case FLOAT: ((float*)vec->array)[++vec->index] = *(float*)elem;
      break;
    case DOUBLE: ((double*)vec->array)[++vec->index] = *(double*)elem;
      break;
    case BOOL: ((bool*)vec->array)[++vec->index] = *(bool*)elem;
      break;
    default:
      // This should never happen, but its an error either way.
      exit(1);
  }
}

void vectorResize(Vector* vec) {
  assert(vec != NULL);
  vec->bufSize = (int)(1.25*vec->bufSize*vec->bytesPerElement);
  realloc(vec->array, vec->bufSize);
}

void vectorTrim(Vector* vec) {
  assert(vec != NULL);
  vec->bufSize = vec->index+1;
  void* arr = malloc(vec->bufSize*vec->bytesPerElement);
  memcpy(arr, vec->array, vec->bufSize*vec->bytesPerElement);;
  if(vec->callbackFree) {
    vec->callbackFree(vec->array);
  } else {
    free(vec->array);
  }
  vec->array = arr;
}
