#ifndef SIMPLEX_H
#define SIMPLEX_H

extern const size_t K;

void transform(const double th[K], double y[K-1]);
void inverse_transform(const double y[K-1], double th[K]);

#endif
