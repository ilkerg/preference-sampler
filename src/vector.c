#include "vector.h"

#include <stdlib.h>

vector_double * vector_double_new(unsigned int size) {
    vector_double *v = malloc(sizeof(vector_double));
    v->data = malloc(size * sizeof(double));
    v->size=size;
    return v;
}

void vector_double_free(vector_double *v) {
    free(v);
}

