
#pragma once 

#include <core.h>

static vec3s const world_forward = { 0.0f, 0.0f, -1.0f };
static vec3s const world_right = { 1.0f, 0.0f, 0.0f };
static vec3s const world_up = { 0.0f, -1.0f, 0.0f };

typedef struct ray 
{
    vec3s origin;
    vec3s direction;
} ray;

void gard_distribution_rays(ray base_ray, float roughness, int num_rays, float frame_rotation, ray *out_rays);

typedef struct plane 
{
    vec3s point;
    vec3s normal;
} plane;

bool ray_plane_intersection(ray r, plane p, float *out_distance);

typedef struct sphere 
{
    vec3s position;
    float radius;
} sphere;

bool ray_sphere_intersection(ray r, vec3s center, float radius, float *out_distance);

typedef enum axis 
{
    axis_x,
    axis_y,
    axis_z
} axis;

typedef struct aabb 
{
    vec3s min, max;
} aabb;

bool ray_aabb_intersection(aabb box, ray r, float *out_distance);

/** returns an empty aabb */
aabb aabb_empty();

/** add a point into the given aabb, growing it */
void aabb_add_point(aabb *box, vec3s point);

/** merge two aabb's into one */
aabb aabb_merge(aabb a, aabb b);

/** return the longest axis of the aabb */
axis aabb_longest_axis(aabb box);


// Convert an sRGB color (in [0,1]) to linear space.
vec3s srgb_to_linear(vec3s c);

// Convert a linear color (in [0,1]) to sRGB space.
vec3s linear_to_srgb(vec3s c);