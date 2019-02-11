#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>

typedef struct matrix {
    size_t nrow;
    size_t ncol;
    double *data;
} matrix;

matrix* matrix_new(size_t nrow, size_t ncol);
void matrix_free(matrix *m);


#endif
