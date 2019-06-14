#ifndef SET_COUNTER_H
#define SET_COUNTER_H

#include <stdbool.h>
#include <stddef.h>

#define HASH_TABLE_SIZE 1024

extern const size_t L;

struct set_counter
{
    bool occupied[HASH_TABLE_SIZE];
    size_t *keys;
    size_t *values;
    size_t size;
};


struct set_counter * set_counter_alloc();
void set_counter_free(struct set_counter *c);
void set_counter_add(struct set_counter *c, const size_t key[L]);
void set_counter_keys(const struct set_counter *c, size_t (*keys)[L]);
void set_counter_values(const struct set_counter *c, size_t *values);

#endif

