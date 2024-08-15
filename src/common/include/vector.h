#ifndef VECTOR_H
#define VECTOR_H

enum DataType {INT, LONG, FLOAT, DOUBLE, BOOL, CHAR};
typedef void(*CallbackFree)(void *);
typedef struct Vector {
  enum DataType dataType;
  int bytesPerElement;
  int bufSize;
  int index;
  void* array;
  CallbackFree callbackFree;
} Vector;

Vector* vectorCreate(int bytesPerElement, int initialCapacity, CallbackFree cbFree, enum DataType dt);
Vector* vectorFromArray(int bytesPerElement, int newCapacity, CallbackFree cbFree, void* initArray);
void freeVector(Vector* vec);
void vectorAppend(Vector* vec, void* elem);
void vectorResize(Vector* vec);

#endif //VECTOR_H
