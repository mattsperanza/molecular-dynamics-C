// Author(s): Matthew Speranza
#ifndef TENSOR_H
#define TENSOR_H

#include "../system/system.h"

void ewaldSource(REAL* src, REAL* r, REAL beta, int tensorOrder);
void tholeSource(REAL* src, REAL* r, REAL ai, REAL aj);
void generateTensor(REAL* tensor, REAL* r, REAL* src, int order);
void multipoleInteraction(System* system, int i, int j, REAL* r, REAL* elecMask);

#endif //TENSOR_H
