
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

float snap_float(float x, float precision)
{
    return roundf(x / precision) * precision;
}

uint64_t ceil_to_interval(uint64_t interval, uint64_t x)
{
    return interval * ceil((double)x / (double)interval);
}


void array_create(array *array, uint64_t type_size, uint64_t page_size)
{
    array->type_size = type_size;
    array->page_size = page_size;

    array->capacity = 0;
    array->size = 0;

    array->data = NULL;
}

void array_destroy(array *array)
{
    array->type_size = 0;
    array->page_size = 0;

    array->capacity = 0;
    array->size = 0;

    if(array->data)
    {
        free(array->data);
    }
}

void array_reserve(array *array, uint64_t preferred_capacity)
{
    uint64_t new_capacity = ceil_to_interval(array->page_size, preferred_capacity);
    if(array->capacity >= new_capacity)
        return;

    void *new_data = malloc(new_capacity * array->type_size);

    if(array->data)
        memcpy(new_data, array->data, array->size * array->type_size);

    free(array->data);
    array->data = new_data;
    array->capacity = new_capacity;
}

void array_shrink(array *array)
{
    if(!array->data)
        return;

    uint64_t new_capacity = ceil_to_interval(array->page_size, array->size);
    if(array->capacity == new_capacity)
        return;

    void *new_data = malloc(new_capacity * array->type_size);
    memcpy(new_data, array->data, array->size * array->type_size);
    free(array->data);
    array->data = new_data;
    array->capacity = new_capacity;
}


void array_clear(array *array)
{
    array->size = 0;
}

void* array_at(array *array, uint64_t index)
{
    return &((uint8_t*)array->data)[array->type_size * index];
}

void* array_set(array *array, uint64_t index, void const *element)
{
    void *p = array_at(array, index);
    memcpy(p, element, array->type_size);
    return p;
}

void array_swap(array *array, uint64_t a, uint64_t b, void *swap_buffer)
{
    void *pa = array_at(array, a);
    void *pb = array_at(array, b);
    memcpy(swap_buffer, pa, array->type_size);
    memcpy(pa, pb, array->type_size);
    memcpy(pb, swap_buffer, array->type_size);
}

void* array_push(array *array, void const *element)
{
    array_reserve(array, array->size + 1);
    void* new_element = array_at(array, array->size);
    memcpy(new_element, element, array->type_size);
    array->size++;
    return new_element;
}

void* array_pop(array *array)
{
    if(array->size == 0)
        return NULL;

    array->size--;
    return array_at(array, array->size);
}

void* array_first(array *array)
{
    return array_at(array, 0);
}

void* array_last(array *array)
{
    return array_at(array, array->size - 1);
}

void* array_remove_fast(array *array, uint64_t index, void *swap_buffer)
{
    if(array->size == 0)
        return NULL;
    
    if(index != array->size - 1 && array->size > 1)
        array_swap(array, index, array->size - 1, swap_buffer);
    
    return array_pop(array);
}

void* array_remove_ordered(array *array, uint64_t index, void *swap_buffer)
{
    if(array->size == 0)
        return NULL;

    if(array->size == 1 || index >= array->size - 1)
        return array_pop(array);

    uint8_t *p = array_at(array, index);
    uint8_t *lp = array_last(array);
    memcpy(swap_buffer, p, array->type_size);
    memmove(p, p + array->type_size, (array->size - index - 1) * array->type_size);
    memcpy(lp, swap_buffer, array->type_size);
    array->size--;
    return lp;
}

void* array_insert(array *array, uint64_t index, void const *element)
{
    if(array->size == 0 || index >= array->size)
        return array_push(array, element);

    array_reserve(array, array->size + 1);

    uint8_t *t = array_at(array, index);
    memmove(t + array->type_size, t, (array->size - index) * array->type_size);
    memcpy(t, element, array->type_size);
    array->size++;
    return t;
}

uint64_t array_index_of(array *array, void const *element)
{
    uint8_t const* element_data = element;
    uint8_t* array_data = array->data;
    return (element_data - array_data) / array->type_size;
}

int64_t array_bsearch(array* array, void const *key, compare_func cmp_func, int64_t *ordered_index)
{

    if(array->size == 0)
    {
        return -1;
    }

    int64_t low = 0;
    int64_t high = array->size - 1;
    int64_t mid = 0;
    int32_t res = 0;
    while(low <= high)
    {
        mid = low + (high - low) / 2;
        void *elem = array_at(array, mid);
        res = cmp_func(key, elem);
        
        if(res == 0)
            return mid;
        
        if (res > 0)
            low = mid + 1;
        else
            high = mid - 1;
    }

    if(ordered_index)
    {
        if(res > 0)
            *ordered_index = mid + 1;
        else 
            *ordered_index = mid;
    }

    return -1;
}

void array_qsort(array* array, compare_func cmp_func)
{
    qsort(array->data, array->size, array->type_size, cmp_func);
}


void* array_insert_ordered(array *array, bool unique, void const *element, compare_func cmp_func)
{
    if(array->size == 0)
    {
        return array_push(array, element);
    }

    uint64_t new_index = -1;
    uint64_t curr = array_bsearch(array, element, cmp_func, &new_index);

    if(curr != -1)
    {
        if(unique)
            return array_at(array, curr);

        return array_insert(array, curr, element);
    }

    return array_insert(array, new_index, element);
}


uint64_t set_segment_index(struct set *set, uint64_t hash)
{
    uint64_t segment_size = UINT64_MAX / set->num_segments;
    return hash / segment_size;
}

void set_create(struct set *set, uint64_t num_segments, uint64_t key_size, hash_func hash_func, compare_func cmp_func)
{
    memset(set, 0, sizeof(struct set));
    set->key_size = key_size;
    set->cmp_func = cmp_func;
    set->hash_func = hash_func;

    set->swap_buffer = malloc(key_size);

    set->num_segments = num_segments;

    set->segments = malloc(sizeof(struct array) * num_segments);
    for(uint64_t i = 0; i < num_segments; i++)
    {
        array_create(&set->segments[i], key_size, 8);
    }
}


void set_destroy(struct set *set)
{
    free(set->swap_buffer);
    for(uint64_t i = 0; i < set->num_segments; i++)
    {
        array_destroy(&set->segments[i]);
    }

    free(set->segments);
}


void set_compress(struct set *set)
{
    for(uint64_t i = 0; i < set->num_segments; i++)
    {
        array_shrink(&set->segments[i]);
    }  
}

bool set_contains(struct set *set, void *key)
{
    uint64_t hash = set->hash_func(set->key_size, key);
    uint64_t segment_index = set_segment_index(set, hash);

    struct array *segment = &set->segments[segment_index];

    if(segment->size == 0)
        return false;

    if(segment->size == 1)
        return set->cmp_func(key, segment->data) == 0;

    int64_t idx = array_bsearch(segment, key, set->cmp_func, NULL);

    return idx != -1;
}


void set_insert(struct set *set, void *key)
{
    uint64_t hash = set->hash_func(set->key_size, key);
    uint64_t segment_index = set_segment_index(set, hash);

    struct array *segment = &set->segments[segment_index];
    array_insert_ordered(segment, true, key, set->cmp_func);
}


void set_remove(struct set *set, void *key)
{
    uint64_t hash = set->hash_func(set->key_size, key);
    uint64_t segment_index = set_segment_index(set, hash);

    struct array *segment = &set->segments[segment_index];

    uint64_t i = array_bsearch(segment, key, set->cmp_func, NULL);

    array_remove_ordered(segment, i,  set->swap_buffer);
}

void set_clear(struct set *set)
{
    for(uint64_t i = 0; i < set->num_segments; i++)
    {
        struct array *segment = &set->segments[i];
        array_clear(segment);
    }
}


uint64_t set_count(struct set *set)
{
    uint64_t c = 0;
    for(uint64_t i = 0; i < set->num_segments; i++)
    {
        struct array *segment = &set->segments[i];
        c += segment->size;
    }

    return c;
}







uint64_t map_segment_index(struct map *map, uint64_t hash)
{
    uint64_t segment_size = UINT64_MAX / map->num_segments;
    return hash / segment_size;
}

void map_create(struct map *map, uint64_t num_segments, uint64_t key_size, uint64_t value_size, hash_func hash_func, compare_func cmp_func)
{
    memset(map, 0, sizeof(struct set));
    map->key_size = key_size;
    map->value_size = value_size;
    map->cmp_func = cmp_func;
    map->hash_func = hash_func;

    map->swap_buffer = malloc(key_size + value_size);

    map->num_segments = num_segments;

    map->segments = malloc(sizeof(struct array) * num_segments);
    for(uint64_t i = 0; i < num_segments; i++)
    {
        array_create(&map->segments[i], key_size +value_size, 8);
    }
}

void map_destroy(struct map *map)
{
    free(map->swap_buffer);
    for(uint64_t i = 0; i < map->num_segments; i++)
    {
        array_destroy(&map->segments[i]);
    }

    free(map->segments);
}

void map_compress(struct map *map)
{
    for(uint64_t i = 0; i < map->num_segments; i++)
    {
        array_shrink(&map->segments[i]);
    }  
}

void *map_at(struct map *map, void *key)
{
    uint64_t hash = map->hash_func(map->key_size, key);
    uint64_t segment_index = map_segment_index(map, hash);

    struct array *segment = &map->segments[segment_index];

    if(segment->size == 0)
        return NULL;

    if(segment->size == 1)
    {
        bool exists = map->cmp_func(key, segment->data) == 0;
        if(exists)
            return (uint8_t*)segment->data + map->key_size;
        else 
            return NULL;
    }

    int64_t idx = array_bsearch(segment, key, map->cmp_func, NULL);
    bool exists = idx != -1;
    if(exists)
        return (uint8_t*)array_at(segment,idx) + map->key_size;
    else 
        return NULL;
}

void* map_set(struct map *map, void *key, void *value)
{
    uint64_t hash = map->hash_func(map->key_size, key);
    uint64_t segment_index = map_segment_index(map, hash);

    struct array *segment = &map->segments[segment_index];


    // no elements, just prepare and push!
    if(segment->size == 0)
    {
        memcpy(map->swap_buffer, key, map->key_size);
        memcpy((uint8_t*)map->swap_buffer + map->key_size, value, map->value_size);
        return (uint8_t*)array_push(segment, map->swap_buffer) + map->key_size;
    }

    // we have exactly one element, check if it matches and if so just update value and return
    if(segment->size == 1)
    {
        bool exists = map->cmp_func(key, segment->data) == 0;
        if(exists)
        {
            memcpy((uint8_t*)segment->data + map->key_size, value, map->value_size);
            return (uint8_t*)segment->data + map->key_size;
        }
    }

    // find index of existing key
    int64_t new_index = -1;
    int64_t idx = array_bsearch(segment, key, map->cmp_func, &new_index);
    if(idx != -1)
    {
        // we found it, just update value and return 
        void* valptr = (uint8_t*)array_at(segment,idx) + map->key_size;
        memcpy(valptr, value, map->value_size);
        return valptr;
    }

    // didnt exist, so prepare and insert!
    memcpy(map->swap_buffer, key, map->key_size);
    memcpy((uint8_t*)map->swap_buffer + map->key_size, value, map->value_size);
    return (uint8_t*)array_insert(segment, new_index, map->swap_buffer) + map->key_size;
}

void map_remove(struct map *map, void *key)
{
    uint64_t hash = map->hash_func(map->key_size, key);
    uint64_t segment_index = map_segment_index(map, hash);

    struct array *segment = &map->segments[segment_index];

    uint64_t i = array_bsearch(segment, key, map->cmp_func, NULL);

    array_remove_ordered(segment, i,  map->swap_buffer);

}

void map_clear(struct map *map)
{
    for(uint64_t i = 0; i < map->num_segments; i++)
    {
        struct array *segment = &map->segments[i];
        array_clear(segment);
    }
}

uint64_t map_count(struct map *map)
{
    uint64_t c = 0;
    for(uint64_t i = 0; i < map->num_segments; i++)
    {
        struct array *segment = &map->segments[i];
        c += segment->size;
    }

    return c;
}

