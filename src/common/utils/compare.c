// Author(s): Matthew Speranza
#include "../include/compare.h"

int compareInt(const void* a, const void* b) {
    return *(int*)a - *(int*)b;
}

int compareReal(const void* a, const void* b){
    REAL c = *(REAL*)a - *(REAL*)b;
    if (c < 0){
        return -1;
    } else if(c > 0){
        return 1;
    }
    return 0;
}
