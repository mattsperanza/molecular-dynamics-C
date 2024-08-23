// Author(s): Matthew Speranza

#include "system.h"

#include <assert.h>

REAL pbc(REAL x, REAL axisLen){
 assert(axisLen > 0);
 while(x > axisLen/2 || x <= -axisLen/2)
  x = x > 0 ? x - axisLen : x + axisLen;
 return x;
}
