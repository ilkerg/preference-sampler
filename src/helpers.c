#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

double
sum(const double *x, const size_t length)
{
    double sum = 0.0;
    for(size_t i = 0; i < length; i++)
        sum += x[i];

    return sum;
}

double
max(const double *x, const size_t length)
{
    double max = GSL_NEGINF;
    for(size_t i = 0; i < length; i++) {
        if (x[i] > max)
            max = x[i];
    }

    return max;
}

void
ones(double *arr, const size_t length)
{
    for(size_t i = 0; i<length; i++)
        arr[i] = 1.0;
}

void
zeros(double *arr, const size_t length)
{
    memset(arr, 0, length * sizeof(double));
}

double
log_sum_exp(const double *x, size_t length)
{
    size_t i;
    double sum = 0;
    double m = max(x, length);

    for (i = 0; i < length; i++) {
        sum += exp(x[i] - m);
    }

    return m + log(sum);
}

int
is_in(size_t e, const size_t *x, size_t length)
{
    for (size_t i=0; i < length; i++) {
        if (x[i] == e)
            return 1;
    }

    return 0;
}

