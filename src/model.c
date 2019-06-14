#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include "model.h"
#include "helpers.h"

double
fullcond(const size_t comp, const double logtheta[K], size_t ngames,
                const size_t games[ngames][L], const size_t game_counts[ngames],
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

        for (size_t l=0; l<L; l++) {
            if (games[i][l]==comp)
                comp_plays = true;
            else if (games[i][l]==K-1)
                K_plays = true;
        }

        if (comp_plays != K_plays) {
            double logth_game[L];
            for (size_t l = 0; l < L; l++) {
                logth_game[l] = logtheta[games[i][l]];
            }

            ll -= game_counts[i] * log_sum_exp(logth_game, L);
        }
    }

    return ll;
}

