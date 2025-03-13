
#include "texture.h"

#include "stb_image.h"

bool texture_load(texture *tex, char const* path)
{

    stbi_ldr_to_hdr_gamma(1.0f);
    stbi_ldr_to_hdr_scale(1.0f);

    stbi_hdr_to_ldr_gamma(1.0f);
    stbi_hdr_to_ldr_scale(1.0f);

    int32_t x, y, c;
    tex->pixels = stbi_loadf(path, &x, &y, &c, 4);
    tex->width = x;
    tex->height = y;
    return tex->pixels != NULL;
}

void texture_destroy(texture *tex)
{
    free(tex->pixels);
}

vec4s texture_sample_pixel(texture *tex, ivec2s pixel)
{
    int64_t index = pixel.x + tex->width * pixel.y;
    if(index < 0)
        index = 0;
    
    if(index >= tex->width * tex->height)
        index = tex->width * tex->height - 1;

    return tex->pixels[index];
}

vec4s texture_sample_uv(texture *tex, vec2s uv)
{
    float x_f = tex->width * uv.x;
    float y_f = tex->height * uv.y;

    float x_min = floorf(x_f);
    float x_max = x_min + 1.0f;
    float x_delta = x_f - x_min;
    
    float y_min = floorf(y_f);
    float y_max = y_min + 1.0f;
    float y_delta = y_f - y_min;

    vec4s tl = texture_sample_pixel(tex, (ivec2s){ x_min, y_min });
    vec4s tr = texture_sample_pixel(tex, (ivec2s){ x_max, y_min });
    vec4s bl = texture_sample_pixel(tex, (ivec2s){ x_min, y_max });
    vec4s br = texture_sample_pixel(tex, (ivec2s){ x_max, y_max });

    vec4s top = glms_vec4_mix(tl, tr, x_delta);
    vec4s bottom = glms_vec4_mix(bl, br, x_delta);

    return glms_vec4_mix(top, bottom, y_delta);
}

vec4s texture_sample_uvw(texture *tex, vec3s uvw)
{
    vec2s const inv_atan = {0.1591f, 0.3193f};
    vec2s uv = (vec2s){ atan2f(uvw.z, uvw.x), asinf(uvw.y) };
    uv = glms_vec2_mul(uv, inv_atan);
    uv = glms_vec2_adds(uv, 0.5f);

    return texture_sample_uv(tex, uv);
}

