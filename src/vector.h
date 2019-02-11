#ifndef VECTOR_H
#define VECTOR_H

#include <stdlib.h>

typedef struct vector {
    size_t size;
    double *data;
} vector;

vector* vector_new(size_t size);
void vector_free(vector *v);

#endif

