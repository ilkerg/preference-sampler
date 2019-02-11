#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <stdlib.h>

void transform(const size_t K, const double th[K], double y[K-1]);
void inverse_transform(const size_t K, const double y[K-1], double th[K]);

#endif
