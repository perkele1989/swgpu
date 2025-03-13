
#pragma once 

#include <core.h>

/** A texture, can be used as a 2D texture, or as a fake cubemap texture using equirectangular projection */
typedef struct texture 
{
    uint32_t width;
    uint32_t height;
    vec4s *pixels;
} texture;

bool texture_load(texture *tex, char const* path);
void texture_destroy(texture *tex);

vec4s texture_sample_pixel(texture *tex, ivec2s pixel);
vec4s texture_sample_uv(texture *tex, vec2s uv);
vec4s texture_sample_uvw(texture *tex, vec3s uvw);
