#include "matrix.h"

#include <stdlib.h>

matrix* matrix_new(size_t nrow, size_t ncol) {
    matrix *m = malloc(sizeof(matrix));
    m->data = malloc(nrow * ncol * sizeof(double));
    m->nrow = nrow;
    m->ncol = ncol;

    return m;
}

void matrix_free(matrix *m) {
    free(m->data);
    free(m);
}

