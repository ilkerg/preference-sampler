#ifndef VECTOR_DOUBLE_H
#define VECTOR_DOUBLE_H

typedef struct vector_double {
    unsigned int size;
    double *data;
} vector_double;

vector_double * vector_double_new(unsigned int size);
void vector_double_free(vector_double *vector);

#endif
