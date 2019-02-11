#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include "helpers.h"

/*
 * S: number of games
 * M: number of players in each game (subset size)
 * K: total number of players
 * x[s]: M-dimensional vector containing player indices of game s
 * y[s]: winner of game s
 */

double loglik(const size_t S, const size_t M, const size_t K,
        const double theta[K], const size_t x[S][M], const size_t y[S]) {
    size_t s, m;
    double ll = 0.;
    double sum_theta;

    for (s = 0; s < S; s++) {
        sum_theta = 0.;
        for (m = 0; m < M; m++) {
            sum_theta += theta[x[s][m]];
        }

        ll += gsl_sf_log(theta[y[s]]) - gsl_sf_log(sum_theta);
    }

    return ll;
}

double fullcond(const size_t S, const size_t M, const size_t K, const size_t comp,
        const double theta[K], const size_t x[S][M], const size_t y[S]) {
    double ll = .0;
    double sum_theta;

    for (size_t s=0; s<S; s++) {
        if (y[s] == comp) {
            ll += gsl_sf_log(theta[comp]);
        }
        else if (y[s] == K-1) {
            ll += gsl_sf_log(theta[K-1]);
        }

        if (is_in(comp, x[s], M) != is_in(K-1, x[s], M)) {
            sum_theta = .0;
            for (size_t m = 0; m < M; m++) {
                sum_theta += theta[x[s][m]];
            }

            ll -= gsl_sf_log(sum_theta);
        }
    }
    return ll;
}
