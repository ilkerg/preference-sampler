#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include "model.h"

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

