#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <immintrin.h>

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
        __m256i idx = _mm256_lddqu_si256((__m256i const *)games[i]);
        __m256i vcomp = _mm256_set1_epi64x(comp);
        __m256i vKm1 = _mm256_set1_epi64x(K-1);

        __m256i xx = _mm256_cmpeq_epi64(idx, vcomp);
        __m256i yy = _mm256_cmpeq_epi64(idx, vKm1);
        __m256i zz = _mm256_or_si256(xx, yy);

        if (_mm256_testz_si256(xx, zz) != _mm256_testc_si256(xx, zz)) {
            __m256d th_game_i = _mm256_i64gather_pd(logtheta, idx, 8);

            ll -= game_counts[i] * log_sum_exp((const double *)&th_game_i, M);
        }
    }

    return ll;
}

