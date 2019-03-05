#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include "model.h"
#include "helpers.h"

double
fullcond(const size_t comp, const double theta[K], size_t ngames,
                const size_t games[ngames][M], const size_t game_counts[ngames],
                const size_t win_counts[K])
{
    double ll = .0;

    double lth[K];
#pragma omp simd
    for (size_t k=0; k<K; k++)
        lth[k] = log(theta[k]);

    /* winners */
    if (win_counts[comp] > 0)
        ll += win_counts[comp] * lth[comp];

    if (win_counts[K-1] > 0)
        ll += win_counts[K-1] * lth[K-1];

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
            double lth_game[M];
            for (size_t m = 0; m < M; m++) {
                lth_game[m] = lth[games[i][m]];
            }

            ll -= game_counts[i] * log_sum_exp(lth_game, M);
        }
    }

    return ll;
}

