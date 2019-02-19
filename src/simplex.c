#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>

#include "simplex.h"

void
transform(const double th[K], double y[K-1])
{
    double sum = 0.;

    size_t Km1 = K-1;
    for (size_t k=0; k<Km1; k++) {
        sum += th[k];
        //s -= th[k];
        y[k] = gsl_sf_log(th[k]) - gsl_sf_log_1plusx(-sum) + gsl_sf_log(Km1-k);
    }
}

void
inverse_transform(const double y[K-1], double th[K])
{
    double z;
    // double s = 1.;
    double sum = 0.;
    size_t Km1 = K-1;
    for (size_t k=0; k<Km1; k++) {
        // z = 1. / (1. + (Km1-k) * gsl_sf_exp(-y[k]));
        z = gsl_sf_exp( -gsl_sf_log_1plusx( (Km1-k)* gsl_sf_exp(-y[k])) );
        // th[k] = s * z;
        th[k] = gsl_sf_exp( gsl_sf_log(1. - sum) + gsl_sf_log(z));
        if(gsl_finite(th[k]) == 0) {
            /*
            fprintf(stderr, "\n");
            fprintf(stderr, "th[%zu] is not finite\n", k);
            fprintf(stderr, "y[%zu] = %.30lf\n", k, y[k]);
            fprintf(stderr, "sum = %.30lf\n", sum);
            */
            th[k] = DBL_MIN;
        }
        sum += th[k];
        // s -= th[k];
    }

    /*
    if (sum > 1.0) {
        fprintf(stderr, "weird stuff happens in numeric computation\n");
        fprintf(stderr, "sum = %.20lf\n", sum);
        double norm=0.;
        for (size_t k=0; k<Km1; k++) {
            fprintf(stderr, "%.20lf\n", y[k]);
            norm += y[k];
        }
        fprintf(stderr, "norm of y = %.20lf\n", norm);
        exit(0);
    }

    assert(sum <= 1.0);
    */

    th[Km1] = (1.0 - sum);
    //assert(fabs(smm - 1.0) < 1e-3);
}

