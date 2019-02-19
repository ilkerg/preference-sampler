#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>

#include <getopt.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "helpers.h"
#include "model.h"
#include "simplex.h"


#define to_string(arr, k) \
    for (size_t i = 0; i < k-1; i++) { \
        printf("%lf,", arr[i]); \
    } \
    printf("%lf", arr[k-1]);


void usage() {
    printf("\n");
    printf("Usage: sampler [options]\n");
    printf("    -h Show this help \n");
    printf("    -v Verbose output \n");
    printf("\n");
    printf("    -k INT (default = 3)        number of dimensions of the simplex\n");
    printf("    -n INT (default = 1000)     number of particles\n");
    printf("    -s INT (default = 100)      number of simulation steps\n");
    printf("    -m INT (default = 2)        subset size\n");
    printf("    -a K-vector (default = ones(K)) prior\n");
    printf("\n");
}

struct opt {
    size_t k;
    size_t n;
    size_t s;
    size_t m;
    int verbose;
    double *alpha;
};

struct sim_context {
    const gsl_rng *r;
    int n;
    int k;
    int s;
    int m;
    int verbose;
    double *alpha;
    double *theta_star;
};


unsigned int move_gibbs(const gsl_rng *r, const size_t S, const size_t M, const size_t K,
        const double th[K], double new_th[K], const size_t x[S][M], const size_t y[S]) {
    double theta_p[K];
    double ll, ll_p;
    double alpha = .0;
    double current_total = .0;
    unsigned int accepted = 0;

    memcpy(theta_p, th, sizeof theta_p);
    for (size_t comp=0; comp<K-1; comp++) {
        /* assert(gsl_fcmp(sum(theta_p, K), 1.0, 1e-10)); */
        ll = fullcond(S, M, K, comp, theta_p, x, y);
        /* determine current sum of all the components
         * except `comp` and the last one
         * because the last component is chosen to be
         * implicitly defined by the others
         */
        current_total = sum(theta_p, K-1) - theta_p[comp];
        /* sample a suitable value for the current component */
/* #pragma omp critical */
        {
        theta_p[comp] = gsl_rng_uniform_pos(r) * (1 - current_total);
        }
        assert(theta_p[comp] > .0);
        theta_p[K-1] = 1 - sum(theta_p, K-1);
        assert(theta_p[K-1] > .0);
        /* assert(gsl_fcmp(sum(theta_p, K), 1.0, 1e-10)); */
        /* theta_p[K-1] -= theta_p[comp] - theta[n][comp]; */
        /* compute full conditional density at theta_p */
        ll_p = fullcond(S, M, K, comp, theta_p, x, y);

/* #pragma omp critical */
        {
        alpha = gsl_rng_uniform_pos(r);
        }
        if (gsl_sf_log(alpha) < ll_p - ll) {
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
    memcpy(new_th, theta_p, sizeof theta_p);
    return accepted;
}

unsigned int move_simplex_transform(const gsl_rng *r, const size_t S, const size_t M, const size_t K,
        const double th[K], double new_th[K], const size_t games[S][M], const size_t winners[S]) {
    unsigned int accepted;

    // forward transform
    double y[K-1];
    double y_prime[K-1];
    double th_p[K];
    transform(K, th, y);

    // move
    gsl_vector *mu = gsl_vector_alloc(K-1);
    gsl_vector *res = gsl_vector_alloc(K-1);
    for(size_t k=0; k<K-1; k++)
        gsl_vector_set(mu, k, y[k]);
    gsl_matrix *sigma = gsl_matrix_alloc(K-1, K-1);
    gsl_matrix_set_identity(sigma);
    //gsl_matrix_scale(sigma, 1);

    gsl_ran_multivariate_gaussian(r, mu, sigma, res);
    memcpy(y_prime, res->data, sizeof y_prime);

    gsl_matrix_free(sigma);
    gsl_vector_free(mu);
    gsl_vector_free(res);

    // inverse transform
    inverse_transform(K, y_prime, th_p);

    // accept-reject
    double ll = loglik(S, M, K, th, games, winners);
    double ll_p = loglik(S, M, K, th_p, games, winners);

    double alpha = gsl_rng_uniform_pos(r);
    if (gsl_sf_log(alpha) < ll_p - ll) {
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


unsigned int move_dirichlet(const gsl_rng *r, const size_t S, const size_t M, const size_t K,
        const double th[K], double new_th[K], const size_t games[S][M], const size_t winners[S]) {
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

    // accept-reject
    double ll = loglik(S, M, K, th, games, winners);
    double ll_p = loglik(S, M, K, th_p, games, winners);

    double alpha_p[K];
    for (size_t k=0; k<K; k++) {
        alpha_p[k] = 100*K*th_p[k];
    }

    double t_xy = gsl_ran_dirichlet_lnpdf(K, alpha, th_p);
    double t_yx = gsl_ran_dirichlet_lnpdf(K, alpha_p, th);

    assert(gsl_finite(ll) == 1);
    assert(gsl_finite(ll_p) == 1);

    double a = gsl_rng_uniform_pos(r);
    if (gsl_sf_log(a) < ll_p - ll + t_yx - t_xy) {
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

unsigned int move_hamiltonian(const gsl_rng *r, const size_t S, const size_t M, const size_t K,
        const double th[K], double new_th[K], const size_t games[S][M], const size_t winners[S]) {

    const size_t L = 20;
    double th_p[K];
    unsigned int accepted;

    double p[K-1];
    double grad[K-1];
    double epsilon;
    double e_initial = .0;
    double e_final = 0.;

    memcpy(th_p, th, K*sizeof(double));

    // sample momentum
    for (size_t k=0; k<K-1; k++) {
        p[k] = gsl_ran_gaussian(r, 1.);
        e_initial += .5 * p[k] * p[k];
    }
    double U = potential(S, M, K, th, games, winners);
    e_initial += U;

    /*
    printf("# theta = ");
    to_string(th_p, K);
    printf("\n");
    printf("# potential = %lf\n", U);
    */

    /* leapfrog (epsilon, L) */
    grad_potential(S, M, K, th, games, winners, grad);
    double amg = fabs(max(grad, K-1));

    epsilon = .05 * sqrt(.2 / amg) / K;
    /* printf("# epsilon = %lf\n", epsilon); */

    for (size_t l=0; l<L; l++) {
        for (size_t k=0; k<K-1; k++) {
            /* printf("# grad[%zu] = %lf\n", k, .5*epsilon*grad[k]); */
            p[k] -= .5*epsilon*grad[k];
            th_p[k] += epsilon*p[k];
        }

        // update implicit coordinate of th_p
        th_p[K-1] = 1. - sum(th_p, K-1);

        /*
        printf("# p [%zu] = ", l);
        to_string(p, K-1);
        printf("\n");
        */

        grad_potential(S, M, K, th_p, games, winners, grad);

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

    // final total energy
    for (size_t k=0; k<K-1; k++)
        e_final += .5 * p[k] * p[k];

    U = potential(S, M, K, th_p, games, winners);
    e_final += U;

    /* printf("# potential = %lf\n", U); */

    /*
     * accept-reject
     */


    /*
    printf("# e_initial = %lf\n", e_initial);
    printf("# e_final = %lf\n", e_final);
    printf("# th_p = ");
    to_string(th_p, K);
    printf("\n");

    printf("# p = ");
    to_string(p, K-1);
    printf("\n");
    */

    assert(gsl_finite(e_initial) == 1);
    assert(gsl_finite(e_final) == 1);
    /* assert(gsl_fcmp(e_initial, e_final, 1e-3) == 0); */

    double a = gsl_rng_uniform_pos(r);
    if (gsl_sf_log(a) < e_final - e_initial) {
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
void resample_move(const gsl_rng *r, const size_t N, const size_t K, const size_t S, const size_t M,
        double theta[N][K], const double w[N], const size_t x[S][M], const size_t y[S]) {

    unsigned int cnt[N];
    size_t accepted = 0;

    gsl_ran_multinomial(r, N, N, w, cnt);

    size_t n_new = 0;
    double (*theta_new)[K] = malloc(N * sizeof *theta_new);
    for (size_t n=0; n<N; n++) {
        if (cnt[n] == 0)
            continue;

/* #pragma omp parallel for reduction(+:accepted) */
        for (size_t i=0; i < cnt[n]; i++) {
            accepted += move_gibbs(r, S, M, K, theta[n], theta_new[n_new+i], x, y);
            /* accepted += move_simplex_transform(r, S, M, K, theta[n], theta_new[n_new+i], x, y); */
            /* accepted += move_dirichlet(r, S, M, K, theta[n], theta_new[n_new+i], x, y); */
            /* accepted += move_hamiltonian(r, S, M, K, theta[n], theta_new[n_new+i], x, y); */
        }
        n_new+=cnt[n];
    }

    printf("# to_move = %zu\n", N);
    printf("# accepted = %zu\n", accepted);
    printf("# acceptance ratio = %lf\n", (double) accepted / N);

    memcpy(theta, theta_new, N*K*sizeof(double));
    free(theta_new);
}

void sample_theta_star(const gsl_rng *r, const size_t K, double theta_star[K]) {
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

void sim(struct sim_context *ctx) {
    /*
     * theta stores current particles
     * which represent current posterior density
     *
     */
    const gsl_rng *r = ctx->r;
    size_t N = ctx->n;
    size_t K = ctx->k;
    size_t S = ctx->s;
    size_t M = ctx->m;
    double *alpha = ctx->alpha;
    double *theta_star = ctx->theta_star;

    double (*theta)[K] = malloc(N * sizeof *theta);
    /* double theta_c[K]; */
    size_t (*x)[M] = malloc(S * sizeof *x);
    size_t *y = malloc(S * sizeof(size_t));

    double w[N];
    double logw[N];

    ones(w, N);
    zeros(logw, N);

    /* general info */
    printf("# generator type: %s\n", gsl_rng_name(r));
    printf("# seed = %lu\n", gsl_rng_default_seed);
    printf("\n");

    /* sample N particles from the prior */
    for (size_t n = 0; n < N; n++) {
        gsl_ran_dirichlet(r, K, alpha, theta[n]);
    }

    for(size_t s = 0; s < S; s++) {
        fprintf(stderr, "s = %zu\r", s);
//#pragma omp parallel for
        for(size_t n=0; n<N; n++) {
            double smm = sum(theta[n], K);
            //fprintf(stderr, "%.30lf\n", smm);
            assert(gsl_finite(smm) == 1);
            assert(gsl_fcmp(smm, 1.0, 1e-10) == 0);
            /*
            double smm = sum(theta[n], K);
            fprintf(stderr, "s=%zu, n=%zu\n", s, n);
            fprintf(stderr, "%.16lf\n", smm);
            assert(fabs(smm - 1.0) < 1e-5);
            */
            /*
            double sm = sum(theta[n], K);
            if(fabs(sm - 1.0) > 1e-3) {
                fprintf(stderr, "\nsum of theta[%zu] at iteration %zu is %lf\n", n, s, sm);
                exit(1);
            }
            */
        }
        /*
         * sample a theta from the current posterior
         */

        /* gsl_ran_choose(r, &theta_c, 1, theta, N, sizeof theta_c); */

        /* size_t theta_sample_idx = gsl_rng_uniform_int(r, N); */

        gsl_ran_discrete_t *g = gsl_ran_discrete_preproc(N, w);
        size_t theta_sample_idx = gsl_ran_discrete(r, g);
        gsl_ran_discrete_free(g);
        /*
        if(s == 0) {
            // pick one of the particles uniformly random
            // as they all have equal weight initially
            gsl_ran_choose(r, &theta_c, 1, theta, N, sizeof theta_c);
        } else {
            // sample one of the particles according to their weights
            gsl_ran_discrete_t *g = gsl_ran_discrete_preproc(N, w);
            size_t j = gsl_ran_discrete(r, g);
            memcpy(theta_c, theta[j], sizeof theta_c);
            gsl_ran_discrete_free(g);
        }
        */

        printf("# sampled theta: ");
        /* to_string(theta_c, K); */
        to_string(theta[theta_sample_idx], K);
        printf("\n");

        /*
         *
         * pick M elmeents from theta_c
         * need to employ a strategy here
         *
         */
        size_t players[M];
        double player_w[M];
        /* gsl_sort_largest_index(players, M, theta_c, 1, K); */
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
            logw[n] += gsl_sf_log(theta_winner) - gsl_sf_log(sum_theta);

            w[n] = gsl_sf_exp(logw[n]);
            assert(gsl_finite(w[n]) == 1);
        }

        double two_logw[N];
        for (size_t n=0; n<N; n++)
            two_logw[n] = 2*logw[n];

        double ess = exp(2*log_sum_exp(logw, N) - log_sum_exp(two_logw, N));

        printf("# ess = %lf\n", ess);

        if (ess < .5*N) {
            printf("# resampling at iteration %zu\n", s);
            resample_move(r, N, K, s+1, M, theta, w, x, y);
            ones(w, N);
            zeros(logw, N);
        }

        printf("\n");
    }

    /* resample at the end  */
    printf("# resampling at iteration %zu\n", S);
    resample_move(r, N, K, S, M, theta, w, x, y);
    /* no need to reset the weights at this point but just to be safe */
    ones(w, N);
    zeros(logw, N);

    for(size_t n = 0; n < N; n++) {
        to_string(theta[n], K);
        printf("\n");
    }

    /*
    double smm = 0.;
#pragma omp parallel for reduction(+:smm)
    for(size_t n=0; n<N; n++) {
        smm += sum(theta[n], K);
    }

    fprintf(stderr, "%lf\n", smm);
    */

    free(theta);
    free(x);
    free(y);
}

int parse_arguments(struct opt *opt, int argc, char *argv[]) {
    int c;
    size_t i = 0;

    opterr = 0;

    while( (c = getopt(argc, argv, "vhk:n:s:m:a:")) != -1) {
        switch(c) {
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 'v':
                opt->verbose = 1;
                break;
            case 'k':
                opt->k = atoi(optarg);
                if(opt->k == 0) return 3;
                break;
            case 'n':
                opt->n = atoi(optarg);
                if(opt->n == 0) return 3;
                break;
            case 's':
                opt->s = atoi(optarg);
                if(opt->s == 0) return 3;
                break;
            case 'm':
                opt->m = atoi(optarg);
                if(opt->m == 0) return 3;
                break;
            case 'a':
                opt->alpha = malloc(opt->k * sizeof(double));
                char *token = strtok(optarg, ",");
                while (token != NULL) {
                    opt->alpha[i++] = atof(token);
                    token = strtok(NULL, ",");
                }
                break;
            default:
                return 2;
        }
    }

    if(optind != argc)
        return 1;

    return 0;
}


int main(int argc, char *argv[]) {
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_set_error_handler_off();

    /* parse command line */
    struct opt opt = {
        .k = 3,
        .n = 1000,
        .s = 10,
        .m = 2,
        .verbose = 0,
        .alpha = NULL
    };

    if (parse_arguments(&opt, argc, argv) != 0) {
        usage();
        exit(EXIT_FAILURE);
    }

    if (opt.alpha == NULL) {
        opt.alpha = malloc(opt.k * sizeof(double));
        ones(opt.alpha, opt.k);
    }

    double *theta_star = malloc(opt.k * sizeof(double));
    sample_theta_star(r, opt.k, theta_star);

    printf("# k = %zu\n", opt.k);
    printf("# n = %zu\n", opt.n);
    printf("# s = %zu\n", opt.s);
    printf("# m = %zu\n", opt.m);
    printf("# v = %d\n", opt.verbose);
    printf("# alpha = ");
    to_string(opt.alpha, opt.k);
    printf("\n");
    printf("# theta_star = ");
    to_string(theta_star, opt.k);
    printf("\n");

    struct sim_context ctx;
    ctx.r = r;
    ctx.n = opt.n;
    ctx.k = opt.k;
    ctx.m = opt.m;
    ctx.s = opt.s;
    ctx.alpha = opt.alpha;
    ctx.theta_star = theta_star;

    sim(&ctx);

    gsl_rng_free(r);

    exit(EXIT_SUCCESS);
}

