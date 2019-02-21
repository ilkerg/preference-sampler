#include <stdlib.h>
#include <string.h>

#include "set_counter.h"


/* python hash implementation for tuples */
static size_t
python_tuple_hash(size_t tuple[M]) {
    size_t mult = 1000003UL;
    size_t x = 0x345678UL;
    ssize_t y;

    ssize_t len = M;

    while (--len >=0) {
        y = *tuple++;
        x = (x ^ y) * mult;
        mult += (size_t) (82520UL + len + len);
    }

    x += 97531UL;
    return x;
}

static size_t
hash(const size_t key[M])
{
    size_t s = 0;
    for (size_t k=0; k<M; k++)
        s += key[k];

    return s;
}

static bool
key_equal(const size_t key1[M], const size_t key2[M])
{
    for (size_t i=0; i<M; i++) {
        if (key1[i] != key2[i])
            return false;
    }
    return true;
}

static void insertion_sort(size_t key[M]) {
    for (size_t i=1; i<M; i++) {
        size_t j = i;
        while (j>0 && key[j-1] > key[j]) {
            size_t tmp = key[j-1];
            key[j-1] = key[j];
            key[j] = tmp;
            j--;
        }
    }
}

struct set_counter *
set_counter_alloc()
{
    struct set_counter *c = malloc(sizeof(struct set_counter));
    memset(c->occupied, 0, sizeof c->occupied);
    c->keys = malloc(M * HASH_TABLE_SIZE * sizeof(size_t));
    c->values = malloc(HASH_TABLE_SIZE * sizeof(size_t));
    c->size = 0;
    return c;
}

void
set_counter_free(struct set_counter *c)
{
    free(c->keys);
    free(c->values);
    free(c);
}

void
set_counter_add(struct set_counter *c, size_t key[M])
{
    /* sort the key */
    insertion_sort(key);

    size_t idx = python_tuple_hash(key) % HASH_TABLE_SIZE;

    size_t *kk = c->keys + M*idx;
    while (c->occupied[idx]) {
        if (key_equal(kk, key)) {
            /* key found. increment and return */
            c->values[idx]++;
            return;
        }

        /* collision. check next slot */
        idx++;
        idx %= HASH_TABLE_SIZE;
    }

    /* key not found. add */
    memcpy(kk, key, M*sizeof(size_t));
    c->values[idx] = 1;
    c->occupied[idx] = true;
    c->size++;
}

void set_counter_keys(const struct set_counter *c, size_t (*keys)[M]) {
    size_t idx=0;
    for (size_t i=0; i<HASH_TABLE_SIZE; i++) {
        if (c->occupied[i]) {
            memcpy(keys[idx++], c->keys + i*M, M*sizeof(size_t));
        }
    }
}

void set_counter_values(const struct set_counter *c, size_t *values) {
    for (size_t i=0; i<HASH_TABLE_SIZE; i++) {
        if (c->occupied[i])
            *values++ = c->values[i];
    }
}


