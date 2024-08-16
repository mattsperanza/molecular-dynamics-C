#ifndef VECTOR_H
#define VECTOR_H
#include <stdbool.h>

/**
 * Vector implementation that handles all memory allocations, reallocations, and frees.
 * All the user needs to do is pass in a method for freeing complex datatypes via cbFree,
 * and call vectorCreate/vectorFree.
 */
enum DataType {INT, LONG, FLOAT, DOUBLE, BOOL, CHAR};
typedef void(*CallbackFree)(void *);
typedef struct Vector {
  enum DataType dataType;
  int bytesPerElement;
  int bufSize; // current size capacity
  int size; // aka next index
  void* array;
  CallbackFree callbackFree;
} Vector;

Vector* vectorCreate(int bytesPerElement, int initialCapacity, CallbackFree cbFree, enum DataType dt);
Vector* vectorFromArray(int bytesPerElement, int newCapacity, CallbackFree cbFree, void* initArray);
Vector* vectorCopy(Vector* vec);
void vectorFree(Vector* vec);
void vectorAppend(Vector* vec, void* elem);
void vectorTrim(Vector* vec);
void vectorClear(Vector* vec);
void vectorResize(Vector* vec);

/////////////////////////////////////////// TESTS

void vectorPrint(Vector* vec, int size);
void vectorTest(bool verbose);

#endif //VECTOR_H
