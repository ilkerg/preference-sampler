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

#include "helpers.h"
#include "model.h"


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


void move(const gsl_rng *r, size_t S, size_t M, size_t K, double th[K], double new_th[K], size_t x[S][M], size_t y[S]) {
    double theta_p[K];
    double ll, ll_p;
    double alpha = .0;
    double current_total = .0;

    memcpy(theta_p, th, sizeof theta_p);
    for (size_t comp=0; comp<K-1; comp++) {
        assert(fabs(sum(theta_p, K) - 1.0) < 1e-10);
        ll = fullcond(S, M, K, comp, theta_p, x, y);
        /* determine current sum of all the components
         * except `comp` and the last one
         * because the last component is chosen to be
         * implicitly defined by the others
         */
        current_total = sum(theta_p, K-1) - theta_p[comp];
        /* sample a suitable value for the current component */
        theta_p[comp] = gsl_rng_uniform_pos(r) * (1 - current_total);
        assert(theta_p[comp] > .0);
        theta_p[K-1] = 1 - sum(theta_p, K-1);
        assert(theta_p[K-1] > .0);
        assert(fabs(sum(theta_p, K) - 1.0) < 1e-10);
        /* theta_p[K-1] -= theta_p[comp] - theta[n][comp]; */
        /* compute full conditional density at theta_p */
        ll_p = fullcond(S, M, K, comp, theta_p, x, y);

        alpha = gsl_rng_uniform_pos(r);
        if (gsl_sf_log(alpha) < ll_p - ll) {
            /* accept */
            /* ll = ll_p; */
        } else {
            /* reject */
            /* reset the proposed component back to its original value */
            theta_p[comp] = th[comp];
            theta_p[K-1] = 1 - sum(theta_p, K-1);
            assert(theta_p[K-1] >= 0);
        }
        memcpy(new_th, theta_p, sizeof theta_p);
    }
}

void resample_move_gibbs(const gsl_rng *r, size_t N, int K, int S, int M, double theta[N][K], const double w[N], size_t x[S][M], size_t y[S]) {
    unsigned int cnt[N];

    gsl_ran_multinomial(r, N, N, w, cnt);

    size_t n_new = 0;
    double (*theta_new)[K] = malloc(N * sizeof *theta_new);
    for (size_t n=0; n<N; n++) {
        assert(fabs(sum(theta[n], K) - 1.0) < 1e-10);
        if (cnt[n] == 0)
            continue;

        if (cnt[n] == 1) {
            memcpy(theta_new[n_new++], theta[n], sizeof theta[n]);
        }
        else {
            for (size_t i=0; i < cnt[n]; i++) {
                move(r, S, M, K, theta[n], theta_new[n_new++], x, y);
            }
        }
    }
    memcpy(theta, theta_new, N*K*sizeof(double));
    free(theta_new);
}

void sample_theta_star(const gsl_rng *r, size_t K, double *theta_star) {
    double a[K];
    /* ones(a, K); */
    for (size_t k=0; k<K; k++) {
        a[k] = 1 / (double) K;
    }

    gsl_ran_dirichlet(r, K, a, theta_star);
}

/*
 * alpha: prior for theta
 * theta_star: true preference vector
 *
 */
void sim(struct sim_context *ctx) {
    /*
     * theta stores current particles
     * which represent current posterior density
     *
     * gsl_matrix uses row-major storage
     * each row is an element of the K-simplex
     */
    const gsl_rng *r = ctx->r;
    int N = ctx->n;
    int K = ctx->k;
    int S = ctx->s;
    int M = ctx->m;
    double *alpha = ctx->alpha;
    double *theta_star = ctx->theta_star;

    double (*theta)[K] = malloc(N * sizeof *theta);
    double theta_c[K];
    double theta_mean[K];
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
        for(size_t n=0; n<N; n++) {
            assert(fabs(sum(theta[n], K) - 1.0) < 1e-10);
        }
        /*
         * sample a theta from the current posterior
         */
        gsl_ran_choose(r, &theta_c, 1, theta, N, sizeof theta_c);
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
        to_string(theta_c, K);
        printf("\n");

        /*
         *
         * pick M elmeents from theta_c
         * need to employ a strategy here
         *
         */
        size_t players[M];
        double player_w[M];
        gsl_sort_largest_index(players, M, theta_c, 1, K);

        memcpy(x[s], players, sizeof players);

        printf("# %d largest elements: ", M);
        for (size_t m=0; m<M-1; m++) {
            printf("%zu,", players[m]);
            player_w[m] = theta_star[players[m]];
        }
        printf("%zu\n", players[M-1]);
        player_w[M-1] = theta_star[players[M-1]];

        /* determine outcome using theta_star */
        gsl_ran_discrete_t *g = gsl_ran_discrete_preproc(M, player_w);
        size_t winner = gsl_ran_discrete(r, g);
        gsl_ran_discrete_free(g);

        y[s] = players[winner];
        printf("# winner: %zu\n", y[s]);

        double theta_winner;

        double w_total = 0;

        zeros(theta_mean, K);

        /* update weights */
        for(size_t n = 0; n < N; n++) {
            theta_winner = theta[n][players[winner]];
            double sum_theta = 0;
            for (size_t m=0; m<M; m++) {
                sum_theta += theta[n][players[m]];
            }

            logw[n] += gsl_sf_log(theta_winner) - gsl_sf_log(sum_theta);
            w[n] = gsl_sf_exp(logw[n]);
            assert(gsl_finite(w[n]) == 1);
            w_total += w[n];

            for (size_t k = 0; k < K; k++) {
                theta_mean[k] += w[n] * theta[n][k];
            }
        }

        for (size_t k = 0; k < K; k++) {
            theta_mean[k] /= w_total;
        }

        printf("# mean theta: ");
        to_string(theta_mean, K);
        printf("\n");

        /* resample_only(r, N, K, theta, w); */
        resample_move_gibbs(r, N, K, s+1, M, theta, w, x, y);
        ones(w, N);
        zeros(logw, N);

        /*
         * compute variance of the log weights
         * and perform resampling if necessary
         */

        /*
        double logw_var = gsl_stats_variance(logw, 1, N);
        printf("# variance of log weights: %lf\n", logw_var);
        if (logw_var > 2) {
            resample_move(r, N, K, s, M, theta, w, x, y);
            ones(w, N);
            zeros(logw, N);
        }
        */

        printf("\n");
    }

    for(size_t n = 0; n < N; n++) {
        to_string(theta[n], K);
        printf("\n");
    }

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


int main (int argc, char *argv[]) {
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

