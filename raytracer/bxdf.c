
#include "bxdf.h"
#include "texture.h"
#include <threads.h>


vec3s fresnel_f0(float ior, vec3s albedo, float metallic)
{
    float ior_medium = 1.0f; // air 
    float reflectance = (ior_medium - ior) / (ior_medium + ior);
    reflectance = reflectance * reflectance;

    vec3s f0 = (vec3s){reflectance, reflectance, reflectance};
    f0 = glms_vec3_mix(f0, albedo, metallic);

    return f0;
}

// Schlick's approximation for Fresnel reflectance.
// Returns: F = F0 + (1 - F0) * (1 - cosTheta)^5
vec3s schlick_fresnel(float cosTheta, vec3s f0)
{
    float factor = powf(1.0f - cosTheta, 5.0f);
    vec3s one = { 1.0f, 1.0f, 1.0f };
    vec3s oneMinusF0 = glms_vec3_sub(one, f0);

    return glms_vec3_add(f0, glms_vec3_scale(oneMinusF0, factor));
}

// Returns a random unit vector in the hemisphere around the given normal.
// (A very basic implementation; in production code you’d want a better sampler.)
vec3s random_in_hemisphere(tracer_thread* thread, vec3s normal)
{
    vec3s randVec;

    do
    {
        randVec.x = (genRand(&thread->rand)) * 2.0f - 1.0f;
        randVec.y = (genRand(&thread->rand)) * 2.0f - 1.0f;
        randVec.z = (genRand(&thread->rand)) * 2.0f - 1.0f;
    } while (glms_vec3_dot(randVec, randVec) >= 1.0f);

    randVec = glms_vec3_normalize(randVec);

    // Ensure the random vector is in the same hemisphere as the normal.
    if(glms_vec3_dot(randVec, normal) < 0.0f)
    {
        randVec = glms_vec3_scale(randVec, -1.0f);
    }

    return randVec;
}

/*
 * Computes the direct lighting at the hit point using a Cook–Torrance BRDF.
 * For demonstration purposes, this function assumes a single directional light
 * (e.g. a sun) with a fixed direction, color, and intensity.
 *
 * The BRDF is composed of a diffuse (Lambertian) term and a specular term computed
 * using a GGX normal distribution, a Smith geometry term, and Schlick's Fresnel.
 *
 * Material parameters (in the surface_proxy) are used as follows:
 * - albedo: base color of the surface.
 * - roughness: controls the microfacet distribution (0 = smooth, 1 = very rough).
 * - metallic: interpolates between a dielectric (non-metal) and a metal.
//  */
// vec3s light_cooktorrance(query *q)
// {
//     surface_proxy *surf = q->hit_surface;

//     // Define a fixed directional light.
//     vec3s lightDir   = glms_vec3_normalize((vec3s){ -1.0f, -1.0f, -1.0f });
//     vec3s lightColor = (vec3s){ 1.0f, 1.0f, 1.0f };
//     float lightIntensity = 10.0f;

//     // Compute common vectors.
//     // The view vector points from the hit point toward the camera.
//     vec3s V = glms_vec3_normalize(glms_vec3_scale(q->ray.direction, -1.0f));
//     vec3s L = lightDir;
//     vec3s H = glms_vec3_normalize(glms_vec3_add(V, L));

//     // Dot products (clamped).
//     float NdotL = fmaxf(glms_vec3_dot(q->hit_normal, L), 0.0f);
//     float NdotV = fmaxf(glms_vec3_dot(q->hit_normal, V), 0.0f);
//     float NdotH = fmaxf(glms_vec3_dot(q->hit_normal, H), 0.0f);
//     float VdotH = fmaxf(glms_vec3_dot(V, H), 0.0f);

//     // Base reflectance at normal incidence.
//     // For dielectrics, F0 is typically around 0.04; for metals, F0 equals the albedo.
//     //vec3s F0 = glms_vec3_mix((vec3s){ 0.04f, 0.04f, 0.04f }, surf->albedo, surf->metallic);
//     vec3s f0 = fresnel_f0(surf->ior, surf->albedo, surf->metallic);
//     vec3s F = schlick_fresnel(VdotH, f0);

//     // GGX Normal Distribution Function (NDF).
//     float alpha = surf->roughness * surf->roughness;  // Remapping roughness.
//     float alpha2 = alpha * alpha;
//     float denom = (NdotH * NdotH * (alpha2 - 1.0f) + 1.0f);
//     float D = alpha2 / (GLM_PI * denom * denom);

//     // Geometry function using Smith's method with Schlick–GGX.
//     float k = (surf->roughness + 1.0f);
//     k = (k * k) / 8.0f;
//     float G_V = NdotV / (NdotV * (1.0f - k) + k);
//     float G_L = NdotL / (NdotL * (1.0f - k) + k);
//     float G = G_V * G_L;

//     // Specular term.
//     // Note: division by a small number is guarded with a tiny offset.
//     vec3s specular = glms_vec3_scale(F, (D * G) / (4.0f * NdotV + 0.001f));

//     // Diffuse term (Lambertian).
//     vec3s one = { 1.0f, 1.0f, 1.0f };
//     vec3s kd = glms_vec3_scale(glms_vec3_sub(one, F), 1.0f - surf->metallic);
//     vec3s diffuse = glms_vec3_scale(glms_vec3_mul(kd, surf->albedo), 1.0f / GLM_PI);

//     // Combine the diffuse and specular components.
//     vec3s BRDF = glms_vec3_add(diffuse, specular);
//     vec3s result = glms_vec3_scale(glms_vec3_mul(lightColor, BRDF), NdotL * lightIntensity);

//     ray rr;
//     rr.direction = lightDir;
//     rr.origin = glms_vec3_add(q->hit_location, glms_vec3_scale(q->hit_normal, 0.001f));

//     // query qq = ray_scene_query(rr, q->scene);
//     // if(qq.hit)
//     // {
//     //     result = glms_vec3_zero();
//     // }

//     // result = glms_vec3_zero();

//     if(q->scene->irradiance_func)
//     {
//         ray irradiance_vector;
//         irradiance_vector.direction = q->hit_normal;
//         irradiance_vector.origin = q->hit_location;
//         float irradiance_factor = 1.0f;
//         result = glms_vec3_add(result, glms_vec3_scale( glms_vec3_mul(q->scene->irradiance_func(irradiance_vector), surf->albedo), (1.0f - surf->metallic) * irradiance_factor) );
//     }

//     return result;
// }







// Convert from local (tangent space) to world space using N as 'Z'.
static vec3s to_world_space(vec3s local, vec3s N, vec3s T, vec3s B)
{
    // local is (x=sinTheta*cosPhi, y=sinTheta*sinPhi, z=cosTheta)
    vec3s world = {
        T.x * local.x + B.x * local.y + N.x * local.z,
        T.y * local.x + B.y * local.y + N.y * local.z,
        T.z * local.x + B.z * local.y + N.z * local.z
    };
    return glms_vec3_normalize(world);
}


// Sample a half-vector H in the local (tangent) coordinate system
// for the GGX distribution. alpha = roughness^2.
static vec3s sample_GGX_HalfVector(tracer_thread *thread, float alpha)
{
    float r1 = genRand(&thread->rand);
    float r2 = genRand(&thread->rand);

    float phi      = 2.0f * GLM_PI * r1;

    // Trowbridge-Reitz (GGX) sampling for half-vector
    // A common form is:
    // cosTheta = sqrt((1 - r2)/(1 + (alpha^2 - 1)*r2));
    float cosTheta = sqrtf((1.0f - r2) / (1.0f + (alpha - 1.0f)*r2 + 1e-6f));
    float sinTheta = sqrtf(fmaxf(0.0f, 1.0f - cosTheta * cosTheta));

    vec3s H_local = {
        sinTheta * cosf(phi),
        sinTheta * sinf(phi),
        cosTheta
    };
    return H_local;
}



vec3s material_light(tracer_thread* thread, query* q, vec3s weight, int32_t depth)
{
    return q->hit_surface->emissive;
}


vec3s material_default(tracer_thread* thread, query* q, vec3s weight, int32_t depth)
{
    surface *hit_surf = q->hit_surface;
    tracer_cfg *cfg = thread->tracer->cfg;

    vec2s uv = q->hit_uv;
    uv.y = 1.0 - uv.y;

    vec3s albedo;
    if(hit_surf->albedo_texture)
    {
        albedo = glms_vec4_copy3( texture_sample_uv(hit_surf->albedo_texture, uv) );

        albedo = srgb_to_linear(albedo);
    }
    else 
    {
        albedo = hit_surf->albedo;
    }

    float roughness;
    if(hit_surf->roughness_texture)
    {
        roughness = texture_sample_uv(hit_surf->roughness_texture, uv).r;
    }
    else 
    {
        roughness = hit_surf->roughness;
    }

    roughness = glm_clamp(roughness, 0.01, 1.0f);

    float metallic;
    if(hit_surf->metallic_texture)
    {
        metallic = texture_sample_uv(hit_surf->metallic_texture, uv).r;
    }
    else 
    {
        metallic = hit_surf->metallic;
    }

    vec3s N, T, B;
    T = glms_vec3_normalize(q->hit_tangent);
    B = glms_vec3_normalize(q->hit_bitangent);
    N = glms_vec3_normalize(q->hit_normal);

    if(hit_surf->normal_texture)
    {
        vec3s normal_map = glms_vec4_copy3(texture_sample_uv(hit_surf->normal_texture, uv));
        normal_map = glms_vec3_subs(glms_vec3_scale(normal_map, 2.0f), 1.0f);
        vec3s mapped_normal;
        mapped_normal.x = T.x * normal_map.x +
                        B.x * normal_map.y +
                        N.x * normal_map.z;

        mapped_normal.y = T.y * normal_map.x +
                        B.y * normal_map.y +
                        N.y * normal_map.z;

        mapped_normal.z = T.z * normal_map.x +
                        B.z * normal_map.y +
                        N.z * normal_map.z;

        mapped_normal = glms_vec3_normalize(mapped_normal);

        N = mapped_normal;
    }


    // vec3s preview_normal = N;
    // //preview_normal = glms_vec3_adds(glms_vec3_scale(preview_normal, 0.5f), 0.5f);
    // preview_normal.x = fmaxf(0.0f, preview_normal.x);
    // preview_normal.y = fmaxf(0.0f, preview_normal.y);
    // preview_normal.z = fmaxf(0.0f, preview_normal.z);
    // return preview_normal;

    vec3s V = glms_vec3_normalize(q->ray.direction);
    vec3s Vneg = glms_vec3_scale(V, -1.0f);

    vec3s f0 = fresnel_f0(hit_surf->ior, albedo, metallic);
    float NdotV = fmaxf(glms_vec3_dot(N, Vneg), 0.0f);

    vec3s base_direction = glms_vec3_normalize(glms_vec3_reflect(V, N));

    vec3s origin = glms_vec3_add(q->hit_location, glms_vec3_scale(N, 0.001f));
    ray base_reflection_ray = { origin, base_direction };
    ray base_diffuse_ray = { origin, N };

    float gard_rotation = ((float)thread->sample_index / (float)cfg->num_samples) * (GLM_PI * 2.0);

    int32_t num_reflected_rays = cfg->min_reflected_rays + ((cfg->max_reflected_rays - cfg->min_reflected_rays) * roughness);
 
    float alpha = roughness * roughness;
    vec3s specular_F = glms_vec3_zero();
    for(int32_t i = 0; i < num_reflected_rays; i++)
    {
        vec3s H_local = sample_GGX_HalfVector(thread, alpha);
        vec3s H = to_world_space(H_local, N, T, B);

        vec3s L = glms_vec3_reflect(Vneg, H);

        float straightener = 1.0 - powf(roughness, 0.95f);
        L = glms_vec3_mix(L, base_direction, straightener);
        L = glms_vec3_normalize(L);

        float NdotL = fmaxf(glms_vec3_dot(N, L), 0.0f);
        if(NdotL <= 0.0f)
        {
            continue;
        }

        float NdotH = fmaxf(glms_vec3_dot(N, H), 0.0f);
        float VdotH = fmaxf(glms_vec3_dot(Vneg, H), 0.0f);

        float denom = (NdotH * NdotH * (alpha - 1.0f) + 1.0f);
        denom = denom * denom; // squared
        float D = alpha / fmaxf(GLM_PI * denom, 0.0001f);

        float k = (roughness + 1.0f);
        k = (k * k) / 8.0f;
        float G1_V = NdotV / fmaxf((NdotV * (1.0f - k) + k), 0.0001f);
        float G1_L = NdotL / fmaxf((NdotL * (1.0f - k) + k), 0.0001f);
        float G = G1_V * G1_L;

        vec3s F = schlick_fresnel(VdotH, f0);
        specular_F = glms_vec3_add(specular_F, F);
        float pdf = (D * NdotH) / fmaxf(4.0f * VdotH, 0.0001f);

        float denomBRDF = 4.0f * NdotV * NdotL;
        vec3s specBRDF = glms_vec3_scale(F, (D * G) / fmaxf(denomBRDF, 0.0001f));

        float ratio = (NdotL) / fmaxf(pdf, 0.0001f);
        vec3s reflection_weight = glms_vec3_scale(specBRDF, ratio);

        vec3s final_weight = glms_vec3_scale(reflection_weight, 1.0f / (float)num_reflected_rays);

        ray out_ray;
        out_ray.origin = glms_vec3_add(q->hit_location, glms_vec3_scale(N, 0.001f));
        out_ray.direction = L;

        array_push(&thread->path_stack, &(path_state){ out_ray, glms_vec3_mul(weight, final_weight), depth - 1 });
    }

    static thread_local ray gard_rays[1024];
    specular_F = glms_vec3_scale(specular_F, 1.0 / num_reflected_rays);
    if(metallic < 0.999f)
    {
        gard_distribution_rays(base_diffuse_ray, 1.0f, cfg->num_diffuse_rays, gard_rotation, gard_rays);
        for(int32_t i = 0; i < cfg->num_diffuse_rays; i++)
        {
            float cos_theta = fmaxf(glms_vec3_dot(N, gard_rays[i].direction), 0.0f);
            vec3s diffuse_weight = glms_vec3_scale(
                glms_vec3_mul(albedo, glms_vec3_sub((vec3s){1.0f, 1.0f, 1.0f}, specular_F)),
                ((cos_theta / GLM_PI) / cfg->num_diffuse_rays) * (1.0f - metallic)
            );

            array_push(&thread->path_stack, &(path_state){gard_rays[i], glms_vec3_mul(weight, diffuse_weight), depth - 1});
        }
    }

    return q->hit_surface->emissive;
}
