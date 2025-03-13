

#pragma once 

#include <stdint.h>
#include <stdbool.h>

/** dynamic array */
typedef struct array 
{
    uint64_t type_size;
    uint64_t page_size;

    uint64_t capacity;
    uint64_t size;
    void *data;
} array;


typedef uint64_t (*hash_func)(uint64_t key_size, const void *key);

/** compare function, returns <0 if a is less than b, 0 if a equals b, >0 if a is greater than b */
typedef int32_t (*compare_func)(const void* a, const void* b);

float snap_float(float x, float precision);

/** snaps x to the closest interval that is higher than x */
uint64_t ceil_to_interval(uint64_t interval, uint64_t x);

/** creates array */
void array_create(array *array, uint64_t type_size, uint64_t page_size);

/** destroys array */
void array_destroy(array *array);

/** reserves at least preferred_capacity number of elements */
void array_reserve(array *array, uint64_t preferred_capacity);

/** shrinks capacity to minimum required */
void array_shrink(array *array);

void array_clear(array *array);

/** returns pointer to element at index */
void* array_at(array *array, uint64_t index);

/** sets the data at index by copy */
void* array_set(array *array, uint64_t index, void const *element);

/** swaps the data at the given indices */
void array_swap(array *array, uint64_t a, uint64_t b, void *swap_buffer);

/** pushes a new element by copy and returns a pointer to it */
void* array_push(array *array, void const *element);

/** removes last element and returns a dangling pointer to it */
void* array_pop(array *array);

/** returns first element */
void* array_first(array *array);

/** returns last element */
void* array_last(array *array);

/** removes element at index without preserving order, returns dangling pointer to removed element */
void* array_remove_fast(array *array, uint64_t index, void *swap_buffer);

/** removes element at index and preserves order, returns dangling pointer to removed element */
void* array_remove_ordered(array *array, uint64_t index, void *swap_buffer);

/** inserts element at the given index by copy and returns pointer to it*/
void* array_insert(array *array, uint64_t index, void const *element);

/** returns the index of the given element (has to point to actual array data) */
uint64_t array_index_of(array *array, void const *element);

/** performs a binary search on a properly sorted array, returns index to element that matches or invalid. if no element is found, lower_bound and upper_bound will be set to either the found lower and upper bounds, or invalid if out of range*/
int64_t array_bsearch(array* array, void const *key, compare_func cmp_func, int64_t *ordered_index);

/** performs a quicksort on this array */
void array_qsort(array* array, compare_func cmp_func);

void* array_insert_ordered(array *array, bool unique, void const *element, compare_func cmp_func);





/** */
typedef struct set 
{
    uint64_t key_size;
    compare_func cmp_func;
    hash_func hash_func;

    void * swap_buffer;

    uint64_t num_segments;
    struct array *segments;
} set;

uint64_t set_segment_index(struct set *set, uint64_t hash);

void set_create(struct set *set, uint64_t num_segments, uint64_t key_size, hash_func hash_func, compare_func cmp_func);

void set_destroy(struct set *set);

void set_compress(struct set *set);

bool set_contains(struct set *set, void *key);

void set_insert(struct set *set, void *key);

void set_remove(struct set *set, void *key);

void set_clear(struct set *set);

uint64_t set_count(struct set *set);



/** */
typedef struct map
{
    uint64_t key_size;
    uint64_t value_size;
    compare_func cmp_func;
    hash_func hash_func;

    void * swap_buffer;

    uint64_t num_segments;
    struct array *segments;
} map;
 
uint64_t map_segment_index(struct map *map, uint64_t hash);

void map_create(struct map *map, uint64_t num_segments, uint64_t key_size, uint64_t value_size, hash_func hash_func, compare_func cmp_func);

void map_destroy(struct map *map);

void map_compress(struct map *map);

void *map_at(struct map *map, void *key);

void *map_set(struct map *map, void *key, void *value);

void map_remove(struct map *map, void *key);

void map_clear(struct map *map);

uint64_t map_count(struct map *map);
