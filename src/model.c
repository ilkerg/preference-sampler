#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include "model.h"
#include "helpers.h"

double
fullcond(const size_t comp, const double logtheta[K], size_t ngames,
                const size_t games[ngames][M], const size_t game_counts[ngames],
                const size_t win_counts[K])
{
    double ll = .0;

    /* winners */
    if (win_counts[comp] > 0)
        ll += win_counts[comp] * logtheta[comp];

    if (win_counts[K-1] > 0)
        ll += win_counts[K-1] * logtheta[K-1];

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
            double logth_game[M];
            for (size_t m = 0; m < M; m++) {
                logth_game[m] = logtheta[games[i][m]];
            }

            ll -= game_counts[i] * log_sum_exp(logth_game, M);
        }
    }

    return ll;
}

