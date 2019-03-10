#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_math.h> // gsl_finite

#ifdef _OPENMP
#include <omp.h>
#endif

const size_t N=10000;
const size_t S=10000;
const size_t K=50;
const size_t M=2;

#include "helpers.h"
#include "model.h"
#include "set_counter.h"

#define to_string(arr, k) \
    for (size_t i = 0; i < k-1; i++) { \
        printf("%.16lf,", arr[i]); \
    } \
    printf("%.16lf", arr[k-1]);

unsigned int
move_gibbs(double *restrict random_numbers, double th[K], size_t ngames,
           const size_t games[ngames][M], const size_t game_counts[ngames],
           const size_t win_counts[K])
{
    double ll, ll_p;
    double alpha;
    double th_comp_old;
    unsigned int accepted;

    for (size_t comp=0; comp<K-1; comp++) {
        assert(gsl_fcmp(sum(th, K), 1.0, 1e-15) == 0);
        ll = fullcond(comp, th, ngames, games, game_counts, win_counts);

        /* sample a suitable value for the current component */
        th_comp_old = th[comp];
        th[comp] = *random_numbers++ * (th[comp] + th[K-1]);

        assert(th[comp] > .0);
        th[K-1] += th_comp_old - th[comp];
        assert(th[K-1] > .0);

        /* compute full conditional density at th_p */
        ll_p = fullcond(comp, th, ngames, games, game_counts, win_counts);

        alpha = *random_numbers++;
        if (log(alpha) < ll_p - ll) {
            /* accept */
            accepted = 1;
        }
        else {
            /* reject */
            /* reset the proposed component back to its original value */
            th[K-1] += th[comp] - th_comp_old;
            th[comp] = th_comp_old;
            assert(th[K-1] >= 0);
            accepted |= 0;
        }
    }

    return accepted;
}


void
resample_move(const gsl_rng *r, double theta[N][K], const double w[N],
              const struct set_counter *games_counter, const size_t wins[K])
{
    unsigned int cnt[N];
    size_t accepted = 0;

    gsl_ran_multinomial(r, N, N, w, cnt);

    /* populate particles */
    double (*theta_new)[K] = malloc(N * sizeof *theta_new);
    size_t n_new = 0;
    for (size_t n=0; n<N; n++) {
        for (size_t i=0; i < cnt[n]; i++) {
            memcpy(theta_new[n_new++], theta[n], sizeof *theta_new);
        }
    }

    /* pre-generate random numbers to avoid thread synchronization */
    double (*random_numbers)[2*(K-1)] = malloc(N * sizeof *random_numbers);
    for (size_t n=0; n<N; n++) {
        for (size_t k=0; k<2*(K-1); k++) {
            random_numbers[n][k] = gsl_rng_uniform_pos(r);
        }
    }


    /* extract game counts from the counter */
    size_t ngames = games_counter->size;
    size_t (*games)[M] = malloc(ngames * sizeof *games);
    size_t *game_counts = malloc(ngames * sizeof *game_counts);
    set_counter_keys(games_counter, games);
    set_counter_values(games_counter, game_counts);

#pragma omp parallel for reduction(+:accepted)
    for (size_t n=0; n<N; n++) {
        accepted += move_gibbs(random_numbers[n], theta_new[n], ngames, games, game_counts, wins);
    }

    printf("# to_move = %zu\n", N);
    printf("# accepted = %zu\n", accepted);
    printf("# acceptance ratio = %lf\n", (double) accepted / N);

    memcpy(theta, theta_new, N*K*sizeof(double));
    free(theta_new);
    free(random_numbers);
    free(games);
    free(game_counts);
}

void
sample_theta_star(const gsl_rng *r, double theta_star[K])
{
    double a[K];
    for (size_t k=0; k<K; k++) {
        a[k] = .1;
    }

    gsl_ran_dirichlet(r, K, a, theta_star);
}

void
read_theta_star(const char *file_name, double theta_star[K])
{
    char buf[80];
    FILE *ts = fopen(file_name, "r");

    if (!ts) {
        fprintf(stderr, "error reading %s\n", file_name);
        exit(EXIT_FAILURE);
    }

    for (size_t k = 0; k<K; k++) {
        if (!fgets(buf, 80, ts)) {
            fprintf(stderr, "error reading %s\n", file_name);
            exit(EXIT_FAILURE);
        }
        theta_star[k] = atof(buf);
    };
    fclose(ts);
}

void
sim(const gsl_rng *r, const double theta_star[K])
{
    double (*theta)[K] = malloc(N * sizeof *theta);
    double *w = malloc(N * sizeof(double));
    double *logw = malloc(N * sizeof(double));

    ones(w, N);
    zeros(logw, N);

    size_t *wins = calloc(K, sizeof *wins);
    struct set_counter *games_counter = set_counter_alloc();

    /* general info */
    printf("# generator type: %s\n", gsl_rng_name(r));
    printf("# seed = %lu\n", gsl_rng_default_seed);
    printf("\n");

    /* sample N particles from the `uniform` prior */
    {
        double alpha[K];
        ones(alpha, K);
        for (size_t n = 0; n < N; n++) {
            gsl_ran_dirichlet(r, K, alpha, theta[n]);
        }
    }

    for(size_t s = 0; s < S; s++) {
        fprintf(stderr, "s = %zu\r", s); /* for progress monitoring */
        printf("# iteration = %zu\n", s);

        /* sample a theta from the current posterior */
        gsl_ran_discrete_t *g = gsl_ran_discrete_preproc(N, w);
        size_t theta_sample_idx = gsl_ran_discrete(r, g);
        gsl_ran_discrete_free(g);

        printf("# sampled theta: ");
        to_string(theta[theta_sample_idx], K);
        printf("\n");

        /* pick M elements from current sample */
        size_t players[M];
        gsl_sort_largest_index(players, M, theta[theta_sample_idx], 1, K);

        set_counter_add(games_counter, players);
        printf("# number of unique subsets so far: %zu\n", games_counter->size);

        double player_w[M];

        printf("# %zu largest elements: ", M);
        for (size_t m=0; m<M-1; m++) {
            printf("%zu,", players[m]);
            player_w[m] = theta_star[players[m]];
        }
        printf("%zu\n", players[M-1]);
        player_w[M-1] = theta_star[players[M-1]];

        printf("# player weights = ");
        to_string(player_w, M);
        printf("\n");

        /* determine outcome using theta_star */
        size_t winner;
        {
            g = gsl_ran_discrete_preproc(M, player_w);
            size_t wn = gsl_ran_discrete(r, g);
            winner = players[wn];
            gsl_ran_discrete_free(g);
        }
        printf("# winner: %zu\n", winner);

        wins[winner]++;

        /* update weights */
#pragma omp parallel for
        for(size_t n = 0; n < N; n++) {
            double theta_winner = theta[n][winner];
            double sum_theta = 0;
            for (size_t m=0; m<M; m++) {
                sum_theta += theta[n][players[m]];
            }
            logw[n] += log(theta_winner) - log(sum_theta);
        }

        /* compute w from logw */
        {
            double lse = log_sum_exp(logw, N);
            for (size_t n=0; n<N; n++) {
                w[n] = exp(logw[n] - lse);
                assert(gsl_finite(w[n]) == 1);
            }
        }

        /* compute ess and perform resampling if necessary */
        {
            double two_logw[N];
            for (size_t n=0; n<N; n++)
                two_logw[n] = 2*logw[n];

            double ess = exp(2*log_sum_exp(logw, N) - log_sum_exp(two_logw, N));

            printf("# ess = %lf\n", ess);

            if (ess < .5*N) {
                printf("# resampling at iteration %zu\n", s);
                resample_move(r, theta, w, games_counter, wins);
                ones(w, N);
                zeros(logw, N);
            }
        }

        printf("\n");
    }
    fprintf(stderr, "\n");

    /* resample at the end  */
    printf("# resampling at iteration %zu\n", S);
    resample_move(r, theta, w, games_counter, wins);
    /* no need to reset the weights at this point but just to be safe... */
    ones(w, N);
    zeros(logw, N);

    for(size_t n = 0; n < N; n++) {
        to_string(theta[n], K);
        printf("\n");
    }

    /* cleanup */
    free(theta);
    free(w);
    free(logw);
    free(wins);
    set_counter_free(games_counter);
}


int
main(int argc, char *argv[])
{
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_set_error_handler_off();

    double *theta_star = malloc(K * sizeof(double));
    /* read theta_star from a file */
    if (argc==2)
        read_theta_star(argv[1], theta_star);
    else
        sample_theta_star(r, theta_star);

    printf("# K = %zu\n", K);
    printf("# N = %zu\n", N);
    printf("# S = %zu\n", S);
    printf("# M = %zu\n", M);
    printf("# theta_star = ");
    to_string(theta_star, K);
    printf("\n");

    // perform simulation
    sim(r, theta_star);

    // cleanup
    free(theta_star);
    gsl_rng_free(r);

    exit(EXIT_SUCCESS);
}

