#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "model.h"
#include "helpers.h"

double
potential(const double theta[K], size_t ngames, const size_t games[ngames][M],
          const size_t winners[ngames])
{
    double V = .0;

    for (size_t k=0; k<K; k++)
        V -= .1*log(theta[k]);

    for (size_t i=0; i<ngames; i++) {
        double sum = .0;
        for (size_t m=0; m<M; m++)
            sum += theta[games[i][m]];

        V -= log(theta[winners[i]]);
        V += log(sum);
    }

    return V;
}

void
grad_potential(const double theta[K], size_t ngames,
               const size_t games[ngames][M], const size_t winners[ngames],
               double grad[K-1])
{
    // initialize grad components
    for (size_t k=0; k<K-1; k++)
        grad[k] = -.1 / theta[k] + .1 / theta[K-1];

    for (size_t i=0; i<ngames; i++) {
        size_t winner = winners[i];
        if (winner == K-1)
            for(size_t k=0; k<K-1; k++)
                grad[k] += 1. / theta[winner];
        else
            grad[winner] -= 1. / theta[winner];

        double sum = 0.;
        for (size_t m=0; m<M; m++) {
            sum += theta[games[i][m]];
        }

        bool K_plays = false;
        for (size_t m=0; m<M; m++) {
            if (games[i][m] == K-1)
                K_plays = true;
        }

        if (K_plays == true) {
            for (size_t k=0; k<K-1; k++) {
                if (is_in(k, games[i], M) == 0)
                    grad[k] -= 1. / sum;
            }
        } else {
            for (size_t m=0; m<M; m++)
                grad[games[i][m]] += 1. / sum;
        }
    }
}

double
loglik(const double theta[K], size_t ngames, const size_t x[ngames][M],
       const size_t y[ngames])
{
    double ll = 0.;
    double sum_theta;

    for (size_t i = 0; i < ngames; i++) {
        sum_theta = 0.;
        for (size_t m = 0; m < M; m++) {
            sum_theta += theta[x[i][m]];
        }

        ll += log(theta[y[i]]) - log(sum_theta);
    }

    return ll;
}

double
fullcond(const size_t comp, const double theta[M], size_t ngames,
                const size_t games[ngames][M], const size_t game_counts[ngames],
                const size_t win_counts[K])
{
    double ll = .0;

    /* winners */
    ll += win_counts[comp] * log(theta[comp]);
    ll += win_counts[K-1] * log(theta[K-1]);

    /* games */

    for (size_t i=0; i<ngames; i++) {
        bool comp_plays = false;
        bool K_plays = false;

        for (size_t m=0; m<M; m++) {
            if (games[i][m]==comp)
                comp_plays = true;
            else if (games[i][m]==K-1)
                K_plays = true;
        }

        if (comp_plays != K_plays) {
            double sum_theta = .0;
            for (size_t m = 0; m < M; m++) {
                sum_theta += theta[games[i][m]];
            }
            ll -= game_counts[i] * log(sum_theta);
        }
    }

    return ll;
}

