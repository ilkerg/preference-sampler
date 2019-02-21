#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "../src/model.h"
#include "../src/set_counter.h"

const size_t K=3;
const size_t M=2;
const size_t S=3;

int main() {
    double tol = 1e-15;

    double theta[K] = { .5, .3, .2 };
    size_t games[S][M] = { {0, 1}, {0, 1}, {0, 2} };
    size_t wins[K] = { 2, 1, 0 };

    struct set_counter *games_counter = set_counter_alloc();

    for (size_t s=0; s<S; s++) {
        set_counter_add(games_counter, games[s]);
    }


    /* extract game counts from the counter */
    size_t ngames = games_counter->size;
    size_t (*unique_games)[M] = malloc(ngames * sizeof *unique_games);
    size_t *game_counts = malloc(ngames * sizeof *game_counts);
    set_counter_keys(games_counter, unique_games);
    set_counter_values(games_counter, game_counts);


    double ll0_expected = 2*log(theta[0]) - 2*log(theta[0] + theta[1]);
    double ll0 = fullcond(0, theta, ngames, unique_games, game_counts, wins);
    assert(fabs(ll0_expected-ll0) < tol);

    double ll1_expected = log(theta[1]) - 2*log(theta[0] + theta[1]) - log(1-theta[1]);
    double ll1 = fullcond(1, theta, ngames, unique_games, game_counts, wins);
    assert(fabs(ll1_expected-ll1) < tol);

    free(unique_games);
    free(game_counts);
    set_counter_free(games_counter);

    printf("ALL GOOD!\n");
    return 0;
}


