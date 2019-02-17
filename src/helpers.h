#ifndef _HELPERS_H
#define _HELPERS_H

double sum(const double *x, const size_t length);
double max(const double *x, const size_t length);
void ones(double *arr, const size_t length);
void zeros(double *arr, const size_t length);
double log_sum_exp(const double x[], const size_t length);
int is_in(const size_t e, const size_t *x, const size_t length);

#endif

