#include <stdio.h>
#include <stdlib.h>

#include "../src/set_counter.h"
#include "../src/int_counter.h"

const size_t M=3;

int
main()
{
    struct set_counter *games = set_counter_alloc();
    struct int_counter *winners = int_counter_alloc();

    size_t g[5][M] = {
        {2, 1, 4},
        {4, 2, 0},
        {0, 1, 2},
        {0, 2, 1},
        {1, 2, 0}
    };

    size_t w[5] = {2, 4, 0, 0, 1};

    for (size_t i=0; i<5; i++) {
        set_counter_add(games, g[i]);
        int_counter_add(winners, w[i]);
    }

    size_t *keys = malloc(games->size * M * sizeof(size_t));
    size_t *values = malloc(games->size * sizeof(size_t));

    set_counter_keys(games, keys);
    set_counter_values(games, values);

    printf("%zu distinct games\n", games->size);

    for (size_t i=0; i<games->size; i++) {
        size_t *kk = keys + i*M;
        for (size_t k=0; k<M; k++) {
            printf("%zu ", kk[k]);
        }
        printf("-> %zu\n", values[i]);
    }

    set_counter_free(games);
    free(keys);
    free(values);

    keys = malloc(winners->size * sizeof(size_t));
    values = malloc(winners->size * sizeof(size_t));

    int_counter_values(winners, values);
    int_counter_keys(winners, keys);

    printf("%zu distinct winners\n", winners->size);

    for (size_t i=0; i<winners->size; i++) {
        printf("%zu -> %zu\n", keys[i], values[i]);
    }

    int_counter_free(winners);
    free(keys);
    free(values);
    return 0;
}

