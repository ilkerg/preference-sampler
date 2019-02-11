#include "vector.h"

#include <stdlib.h>

vector* vector_new(size_t size) {
    vector *v = malloc(sizeof(vector));
    v->data = malloc(size * sizeof(double));
    v->size=size;
    return v;
}

void vector_free(vector *v) {
    free(v->data);
    free(v);
}

