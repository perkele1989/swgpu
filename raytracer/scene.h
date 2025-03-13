
#pragma once 

#include <core.h>
#include <utils.h>

#include "common.h"

typedef struct query query;
typedef struct mesh mesh;
typedef struct tracer_thread tracer_thread;
typedef struct texture texture;

typedef vec3s (*material_func)(tracer_thread* thread, query* q, vec3s weight, int32_t depth);

typedef struct surface 
{
    material_func material_func;

    texture *emissive_texture;
    vec3s emissive;
    texture *albedo_texture;
    vec3s albedo;
    texture *roughness_texture;
    float roughness;
    texture *metallic_texture; 
    float metallic;
    float ior;

    texture *normal_texture;

} surface;

typedef struct entity 
{
    mesh *mesh;
    surface *surface;

    mat4s transform;
    mat4s inverse_transform;
    mat3s normal; 
    mat3s inverse_normal;
} entity;

void entity_set_transform(entity *m, vec3s position, versors rotation);

typedef struct scene 
{
    texture* environment_hdr; 
    float environment_intensity;

    /** type: entity */
    struct array entities; 
} scene;

void scene_create(scene *s);
void scene_destroy(scene *s);

typedef struct query
{
    ray ray;
    scene *scene;

    bool hit;

    float hit_distance;
    vec3s hit_location;
    vec3s hit_normal;
    vec3s hit_tangent;
    vec3s hit_bitangent;
    vec2s hit_uv;

    surface *hit_surface;
} query;

query scene_query(scene *scene, ray ray);


typedef struct fast_query
{
    ray ray;
    scene *scene;

    bool hit;

} fast_query;


fast_query scene_fast_query(scene *scene, ray ray);
