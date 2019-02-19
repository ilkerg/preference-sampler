#ifndef MODEL_H
#define MODEL_H

extern const size_t K;
extern const size_t M;

double loglik(const double theta[M], size_t ngames,
              const size_t games[ngames][M], const size_t y[ngames]);

double fullcond(const size_t comp, const double theta[M],
                size_t ngames, const size_t games[ngames][M],
                const size_t y[ngames]);

double potential(const double theta[K], const size_t ngames,
                 const size_t games[ngames][M], const size_t winners[ngames]);

void grad_potential(const double theta[K], size_t ngames,
                    const size_t games[ngames][M], const size_t winners[ngames],
                    double grad[K-1]);

#endif
