#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>

#ifdef _OPENMP
#include <omp.h>
#endif

const size_t N=10000;
const size_t S=200;
const size_t K=10;
const size_t M=3;

#include "helpers.h"
#include "model.h"
#include "simplex.h"


#define to_string(arr, k) \
    for (size_t i = 0; i < k-1; i++) { \
        printf("%lf,", arr[i]); \
    } \
    printf("%lf", arr[k-1]);


unsigned int
move_gibbs(double random_numbers[2*(K-1)], double th[K], size_t ngames,
           const size_t games[ngames][M], const size_t winners[ngames])
{
    double theta_p[K];
    double ll, ll_p;
    double alpha = .0;
    double current_total = .0;
    unsigned int accepted = 0;

    memcpy(theta_p, th, sizeof theta_p);
    for (size_t comp=0; comp<K-1; comp++) {
        /* assert(gsl_fcmp(sum(theta_p, K), 1.0, 1e-10)); */
        ll = fullcond(comp, theta_p, ngames, games, winners);
        /* determine current sum of all the components
         * except `comp` and the last one
         * because the last component is chosen to be
         * implicitly defined by the others
         */
        current_total = sum(theta_p, K-1) - theta_p[comp];
        /* sample a suitable value for the current component */
        theta_p[comp] = random_numbers[comp] * (1 - current_total);

        assert(theta_p[comp] > .0);
        theta_p[K-1] = 1 - sum(theta_p, K-1);
        assert(theta_p[K-1] > .0);

        /* compute full conditional density at theta_p */
        ll_p = fullcond(comp, theta_p, ngames, games, winners);

        alpha = random_numbers[comp+1];
        if (log(alpha) < ll_p - ll) {
            /* accept */
            /* ll = ll_p; */
            accepted |= 1;
        } else {
            /* reject */
            /* reset the proposed component back to its original value */
            theta_p[comp] = th[comp];
            theta_p[K-1] = 1 - sum(theta_p, K-1);
            assert(theta_p[K-1] >= 0);
            accepted |= 0;
        }
    }
    memcpy(th, theta_p, sizeof theta_p);
    return accepted;
}

unsigned int
move_simplex_transform(const gsl_rng *r, const double th[K], double new_th[K],
                       size_t ngames, const size_t games[ngames][M],
                       const size_t winners[ngames])
{
    unsigned int accepted;

    /* forward transform */
    double y[K-1];
    double y_prime[K-1];
    double th_p[K];
    transform(th, y);

    /* move */
    gsl_vector *mu = gsl_vector_alloc(K-1);
    gsl_vector *res = gsl_vector_alloc(K-1);
    for(size_t k=0; k<K-1; k++)
        gsl_vector_set(mu, k, y[k]);
    gsl_matrix *sigma = gsl_matrix_alloc(K-1, K-1);
    gsl_matrix_set_identity(sigma);

    gsl_ran_multivariate_gaussian(r, mu, sigma, res);
    memcpy(y_prime, res->data, sizeof y_prime);

    gsl_matrix_free(sigma);
    gsl_vector_free(mu);
    gsl_vector_free(res);

    /* inverse transform */
    inverse_transform(y_prime, th_p);

    /* accept-reject */
    double ll = loglik(th, ngames, games, winners);
    double ll_p = loglik(th_p, ngames, games, winners);

    double alpha = gsl_rng_uniform_pos(r);
    if (log(alpha) < ll_p - ll) {
        /* accept */
        memcpy(new_th, th_p, K*sizeof(double));
        accepted = 1;
    } else {
        /* reject */
        memcpy(new_th, th, K*sizeof(double));
        accepted = 0;
    }

    return accepted;
}


unsigned int
move_dirichlet(const gsl_rng *r, const double th[K], double new_th[K],
               size_t ngames, const size_t games[ngames][M],
               const size_t winners[ngames])
{
    double th_p[K];
    unsigned int accepted;

    double alpha[K];
    for (size_t k=0; k<K; k++) {
        alpha[k] = 100*K*th[k];
    }

//#pragma omp critical
    gsl_ran_dirichlet(r, K, alpha, th_p);

    for (size_t k=0; k<K; k++) {
        if (th_p[k] == 0)
            th_p[k] = DBL_MIN;
    }

    /* accept-reject */
    double ll = loglik(th, ngames, games, winners);
    double ll_p = loglik(th_p, ngames, games, winners);

    double alpha_p[K];
    for (size_t k=0; k<K; k++) {
        alpha_p[k] = 100*K*th_p[k];
    }

    double t_xy = gsl_ran_dirichlet_lnpdf(K, alpha, th_p);
    double t_yx = gsl_ran_dirichlet_lnpdf(K, alpha_p, th);

    assert(gsl_finite(ll) == 1);
    assert(gsl_finite(ll_p) == 1);

    double a = gsl_rng_uniform_pos(r);
    if (log(a) < ll_p - ll + t_yx - t_xy) {
        /* accept */
        memcpy(new_th, th_p, K*sizeof(double));
        accepted = 1;
    } else {
        /* reject */
        memcpy(new_th, th, K*sizeof(double));
        accepted = 0;
    }
    return accepted;
}

unsigned int
move_hamiltonian(const gsl_rng *r, const double th[K], double new_th[K],
                 size_t ngames, const size_t games[ngames][M],
                 const size_t winners[ngames])
{
    const size_t L = 20;
    double th_p[K];
    unsigned int accepted;

    double p[K-1];
    double grad[K-1];
    double epsilon;
    double e_initial = .0;
    double e_final = 0.;

    memcpy(th_p, th, K*sizeof(double));

    /* sample momentum */
    for (size_t k=0; k<K-1; k++) {
        p[k] = gsl_ran_gaussian(r, 1.);
        e_initial += .5 * p[k] * p[k];
    }
    double U = potential(th, ngames, games, winners);
    e_initial += U;

    /* leapfrog (epsilon, L) */
    grad_potential(th, ngames, games, winners, grad);
    double amg = fabs(max(grad, K-1));

    epsilon = .05 * sqrt(.2 / amg) / K;

    for (size_t l=0; l<L; l++) {
        for (size_t k=0; k<K-1; k++) {
            p[k] -= .5*epsilon*grad[k];
            th_p[k] += epsilon*p[k];
        }

        /* update implicit coordinate of th_p */
        th_p[K-1] = 1. - sum(th_p, K-1);

        grad_potential(th_p, ngames, games, winners, grad);

        for (size_t k=0; k<K-1; k++) {
            p[k] -= .5*epsilon*grad[k];
        }
    }


    /*
     * if new position is out of the simplex
     * because of numerical issues
     * simply reject
     */
    for (size_t k=0; k<K; k++) {
        if (th_p[k] < 0.) {
            memcpy(new_th, th, K*sizeof(double));
            return 0;
        }
    }

    assert(gsl_fcmp(sum(th_p, K), 1.0, 1e-15) == 0);

    /* final total energy */
    for (size_t k=0; k<K-1; k++)
        e_final += .5 * p[k] * p[k];

    U = potential(th_p, ngames, games, winners);
    e_final += U;

    /* accept-reject */
    assert(gsl_finite(e_initial) == 1);
    assert(gsl_finite(e_final) == 1);

    double a = gsl_rng_uniform_pos(r);
    if (log(a) < e_final - e_initial) {
        /* accept */
        memcpy(new_th, th_p, K*sizeof(double));
        accepted = 1;
    } else {
        /* reject */
        memcpy(new_th, th, K*sizeof(double));
        accepted = 0;
    }
    return accepted;

    return accepted;
}

void
resample_move(const gsl_rng *r, double theta[N][K], const double w[N],
              const size_t ngames, const size_t x[ngames][M],
              const size_t y[ngames])
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

#pragma omp parallel for reduction(+:accepted)
    for (size_t n=0; n<N; n++) {
        accepted += move_gibbs(random_numbers[n], theta_new[n], ngames, x, y);
        /* accepted += move_simplex_transform(r, S, M, K, theta[n], theta_new[n_new+i], x, y); */
        /* accepted += move_dirichlet(r, S, M, K, theta[n], theta_new[n_new+i], x, y); */
        /* accepted += move_hamiltonian(r, S, M, K, theta[n], theta_new[n_new+i], x, y); */
    }

    printf("# to_move = %zu\n", N);
    printf("# accepted = %zu\n", accepted);
    printf("# acceptance ratio = %lf\n", (double) accepted / N);

    memcpy(theta, theta_new, N*K*sizeof(double));
    free(theta_new);
    free(random_numbers);
}

void
sample_theta_star(const gsl_rng *r, double theta_star[K])
{
    double a[K];
    double ts[K];
    /* ones(a, K); */
    for (size_t k=0; k<K; k++) {
        a[k] = 10. / K;
    }

    gsl_ran_dirichlet(r, K, a, ts);
    for (size_t k=0; k<K; k++) {
        assert(ts[k] > 0.);
    }
    gsl_sort_largest(theta_star, K, ts, 1, K);
}

void
sim(const gsl_rng *r, const double theta_star[K])
{
    /*
     * theta stores current particles
     * which represent current posterior density
     */
    double (*theta)[K] = malloc(N * sizeof *theta);
    /* double theta_c[K]; */
    size_t (*x)[M] = malloc(S * sizeof *x);
    size_t *y = malloc(S * sizeof(size_t));

    double *w = malloc(N * sizeof(double));
    double *logw = malloc(N * sizeof(double));

    ones(w, N);
    zeros(logw, N);

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
        fprintf(stderr, "s = %zu\r", s);
//#pragma omp parallel for
        for(size_t n=0; n<N; n++) {
            double smm = sum(theta[n], K);
            assert(gsl_finite(smm) == 1);
            assert(gsl_fcmp(smm, 1.0, 1e-10) == 0);
        }

        /* sample a theta from the current posterior */
        gsl_ran_discrete_t *g = gsl_ran_discrete_preproc(N, w);
        size_t theta_sample_idx = gsl_ran_discrete(r, g);
        gsl_ran_discrete_free(g);

        printf("# sampled theta: ");
        to_string(theta[theta_sample_idx], K);
        printf("\n");

        /* pick M elements from theta_c */
        size_t players[M];
        double player_w[M];
        gsl_sort_largest_index(players, M, theta[theta_sample_idx], 1, K);

        memcpy(x[s], players, sizeof players);

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
        g = gsl_ran_discrete_preproc(M, player_w);
        size_t winner = gsl_ran_discrete(r, g);
        gsl_ran_discrete_free(g);

        y[s] = players[winner];
        printf("# winner: %zu\n", y[s]);

        /* update weights */
//#pragma omp parallel for
        for(size_t n = 0; n < N; n++) {
            double theta_winner = theta[n][players[winner]];
            double sum_theta = 0;
            for (size_t m=0; m<M; m++) {
                sum_theta += theta[n][players[m]];
            }
            logw[n] += log(theta_winner) - log(sum_theta);
        }

        /* compute w from logw */
        for (size_t n=0; n<N; n++) {
            w[n] = exp(logw[n]);
            assert(gsl_finite(w[n]) == 1);
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
                resample_move(r, theta, w, s+1, x, y);
                ones(w, N);
                zeros(logw, N);
            }
        }

        printf("\n");
    }

    /* resample at the end  */
    printf("# resampling at iteration %zu\n", S);
    resample_move(r, theta, w, S, x, y);
    /* no need to reset the weights at this point but just to be safe... */
    ones(w, N);
    zeros(logw, N);

    for(size_t n = 0; n < N; n++) {
        to_string(theta[n], K);
        printf("\n");
    }

    free(theta);
    free(x);
    free(y);
    free(w);
    free(logw);
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

