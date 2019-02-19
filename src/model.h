#ifndef MODEL_H
#define MODEL_H

extern const size_t S;
extern const size_t K;
extern const size_t M;

double loglik(const double theta[M], const size_t x[S][M], const size_t y[S]);

double fullcond(const size_t comp, const double theta[M],
                const size_t x[S][M], const size_t y[S]);

double potential(const double theta[K], const size_t games[S][M],
                 const size_t winners[S]);

void grad_potential(const double theta[K], const size_t games[S][M],
                    const size_t winners[S], double grad[K-1]);

#endif
