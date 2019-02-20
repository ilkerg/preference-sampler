#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "../src/model.h"

const size_t K=3;
const size_t M=2;
const size_t S=1;

int main() {
    double tol = 1e-15;

    double theta[K] = { .5, .3, .2 };
    size_t x[S][M] = { {0, 1} };
    size_t y[S] = { 0 };

    double ll0_expected = log(theta[0]) - log(theta[0] + theta[1]);
    double ll0 = fullcond(0, theta, S, x, y);
    assert(fabs(ll0_expected-ll0) < tol);

    double ll1_expected = -log(theta[0] + theta[1]);
    double ll1 = fullcond(1, theta, S, x, y);
    assert(fabs(ll1_expected-ll1) < tol);
}

