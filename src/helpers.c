#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>

double sum(double x[], size_t length) {
    double sum = 0.0;
    for(size_t i = 0; i < length; i++)
        sum += x[i];

    return sum;
}

void ones(double *arr, size_t length) {
    for(size_t i = 0; i<length; i++)
        arr[i] = 1.0;
}

void zeros(double *arr, size_t length) {
    memset(arr, 0.0, length * sizeof(double));
}

double log_sum_exp(double x[], size_t length) {
    size_t i;
    double sum = 0;
    double max = GSL_NEGINF;

    for (i = 0; i < length; i++) {
        if (x[i] > max)
            max = x[i];
    }

    for (i = 0; i < length; i++) {
        sum += exp(x[i] - max);
    }

    return max + log(sum);
}

int is_in(size_t e, size_t *x, size_t length) {
    for (size_t i=0; i < length; i++) {
        if (x[i] == e)
            return 1;
    }

    return 0;
}
