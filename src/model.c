#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>

#include "helpers.h"

/*
 * S: number of games
 * M: number of players in each game (subset size)
 * K: total number of players
 * x[s]: M-dimensional vector containing player indices of game s
 * y[s]: winner of game s
 */

double loglik(size_t S, size_t M, size_t K, double theta[K], size_t x[S][M], size_t y[S]) {
    size_t s, m;
    double ll = 0;
    double sum_theta;

    for (s = 0; s < S; s++) {
        sum_theta = 0.0;
        for (m = 0; m < M; m++) {
            sum_theta += theta[x[s][m]];
        }

        ll += log(theta[y[s]]) - log(sum_theta);
    }

    return ll;
}

double fullcond(size_t S, size_t M, size_t K, size_t comp, double theta[K], size_t x[S][M], size_t y[S]) {
    double ll = .0;
    double sum_theta;

    for (size_t s=0; s<S; s++) {
        if (y[s] == comp) {
            ll += log(theta[comp]);
        }
        else if (y[s] == K-1) {
            ll += log(theta[K-1]);
        }

        if (is_in(comp, x[s], M) != is_in(K-1, x[s], M)) {
            sum_theta = .0;
            for (size_t m = 0; m < M; m++) {
                sum_theta += theta[x[s][m]];
            }

            ll -= log(sum_theta);
        }
    }
    return ll;
}
