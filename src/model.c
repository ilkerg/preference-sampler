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

double potential(const size_t S, const size_t M, const size_t K,
        const double theta[K], const size_t games[S][M], const size_t winners[S]) {
    double V = .0;

    for (size_t k=0; k<K; k++)
        V -= .1*log(theta[k]);

    for (size_t s=0; s<S; s++) {
        double sum = .0;
        for (size_t m=0; m<M; m++)
            sum += theta[games[s][m]];

        V -= log(theta[winners[s]]);
        V += log(sum);
    }

    return V;
}

void grad_potential(const size_t S, const size_t M, const size_t K,
        const double theta[K], const size_t games[S][M], const size_t winners[S],
        double grad[K-1]) {

    // initialize grad components
    for (size_t k=0; k<K-1; k++)
        grad[k] = -.1 / theta[k] + .1 / theta[K-1];

    for (size_t s=0; s<S; s++) {
        size_t winner = winners[s];
        if (winner == K-1)
            for(size_t k=0; k<K-1; k++)
                grad[k] += 1. / theta[winner];
        else
            grad[winner] -= 1. / theta[winner];

        double sum = 0.;
        int K_plays = 0;
        for (size_t m=0; m<M; m++) {
            sum += theta[games[s][m]];
            if (games[s][m] == K-1)
                K_plays = 1;
        }

        if (K_plays == 1) {
            for (size_t k=0; k<K-1; k++) {
                if (is_in(k, games[s], M) == 0)
                    grad[k] -= 1. / sum;
            }
        } else {
            for (size_t m=0; m<M; m++)
                grad[games[s][m]] += 1. / sum;
        }
    }
}

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

        ll += log(theta[y[s]]) - log(sum_theta);
    }

    return ll;
}

double fullcond(const size_t S, const size_t M, const size_t K, const size_t comp,
        const double theta[K], const size_t x[S][M], const size_t y[S]) {
    double ll = .0;

    for (size_t s=0; s<S; s++) {
        size_t winner = y[s];
        if (winner == comp || winner == K-1) {
            ll += log(theta[winner]);
        }

        unsigned int comp_plays = 0;
        unsigned int K_plays = 0;

        for (size_t m=0; m<M; m++) {
            if (x[s][m]==comp)
                comp_plays = 1;
            else if (x[s][m]==K-1)
                K_plays = 1;
        }

        if (comp_plays != K_plays) {
            double sum_theta = .0;
            for (size_t m = 0; m < M; m++) {
                sum_theta += theta[x[s][m]];
            }
            ll -= log(sum_theta);
        }

        /*
        if (is_in(comp, x[s], M) != is_in(K-1, x[s], M)) {
            sum_theta = .0;
            for (size_t m = 0; m < M; m++) {
                sum_theta += theta[x[s][m]];
            }

            ll -= log(sum_theta);
        }
        */
    }
    return ll;
}
