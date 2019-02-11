
double loglik(const size_t S, const size_t M, const size_t K,
        const double theta[M], const size_t x[S][M], const size_t y[S]);
double fullcond(const size_t S, const size_t M, const size_t K, const size_t comp,
        const double theta[M], const size_t x[S][M], const size_t y[S]);

