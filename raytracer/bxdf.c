
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


/*
 * Computes the reflected ray using a microfacet approach.
 * The perfect mirror reflection is computed first, and then it is perturbed by a random
 * offset that is proportional to the surface roughness. (A rough surface produces a blurrier
 * reflection.) The Fresnel term (using Schlick's approximation) is used as the reflectance weight.
//  */
// int32_t reflection_microfacet(query *q, weighted_ray* out_rays)
// {
//     surface_proxy *surf = q->hit_surface;

//     // Normal and view vectors.
//     vec3s N = glms_vec3_normalize(q->hit_normal);
//     vec3s V = glms_vec3_normalize(q->ray.direction);

//     // Compute perfect mirror reflection direction.
//     vec3s reflectedRay = glms_vec3_normalize(glms_vec3_reflect(V, N));

//     // Perturb the reflection direction by mixing in a random direction from the hemisphere.
//     vec3s randomDir = random_in_hemisphere(N);
//     // Roughness = 0 gives perfect mirror reflection; 1 gives a fully diffused reflection.
//     reflectedRay = glms_vec3_normalize(glms_vec3_mix(reflectedRay, randomDir, surf->roughness));

//     // Compute Fresnel reflectance weight.
//     //vec3s F0 = glms_vec3_mix((vec3s){ 0.04f, 0.04f, 0.04f }, surf->albedo, surf->metallic);
//     vec3s f0 = fresnel_f0(surf->ior, surf->albedo, surf->metallic);
//     float cosTheta = fmaxf(glms_vec3_dot(glms_vec3_scale(V, -1.0f), N), 0.0f);
//     vec3s F = schlick_fresnel(cosTheta, f0);

//     // Create the new reflected ray, offsetting the origin slightly to avoid self-intersection.
//     vec3s origin = glms_vec3_add(q->hit_location, glms_vec3_scale(N, 0.001f));
//     ray newRay = { origin, reflectedRay };

//     out_rays[0].ray = newRay;
//     out_rays[0].weight = F; // Per-channel reflectance.

//     return 1;
// }

// int32_t reflection_microfacet_multi_rand(query *q, weighted_ray* out_rays)
// {
//     surface_proxy *surf = q->hit_surface;

//     // Normal and view vectors.
//     vec3s N = glms_vec3_normalize(q->hit_normal);
//     vec3s V = glms_vec3_normalize(q->ray.direction);

//     // Compute Fresnel reflectance weight.
//     //vec3s F0 = glms_vec3_mix((vec3s){ 0.04f, 0.04f, 0.04f }, surf->albedo, surf->metallic);
//     vec3s f0 = fresnel_f0(surf->ior, surf->albedo, surf->metallic);
//     float cosTheta = fmaxf(glms_vec3_dot(glms_vec3_scale(V, -1.0f), N), 0.0f);
//     vec3s F = schlick_fresnel(cosTheta, f0);

//     // Compute perfect mirror reflection direction.
//     vec3s reflected_ray = glms_vec3_normalize(glms_vec3_reflect(V, N));
//     vec3s origin = glms_vec3_add(q->hit_location, glms_vec3_scale(N, 0.001f));

//     // Perturb the reflection direction by mixing in a random direction from the hemisphere.
//     int32_t num_rays = min_reflected_rays + ((max_reflected_rays-min_reflected_rays) * surf->roughness);
//     for(int32_t i = 0; i < num_rays; i++)
//     {
//         vec3s random_dir = random_in_hemisphere(N);
//         // Roughness = 0 gives perfect mirror reflection; 1 gives a fully diffused reflection.
//         vec3s perturbed_ray = glms_vec3_normalize(glms_vec3_mix(reflected_ray, random_dir, surf->roughness));


//         ray newRay = { origin, perturbed_ray };

//         out_rays[i].ray = newRay;
//         out_rays[i].weight = glms_vec3_scale(F, 1.0f / (float)num_rays);
//     }

//     return num_rays;
// }

// /** same as above, except multireflection using gard distribution */
// int32_t reflection_microfacet_multi_gard(query *q, weighted_ray* out_rays)
// {
//     surface_proxy *surf = q->hit_surface;

//     // Normal and view vectors.
//     vec3s N = glms_vec3_normalize(q->hit_normal);
//     vec3s V = glms_vec3_normalize(q->ray.direction);

//     // Compute Fresnel reflectance weight.
//     vec3s f0 = fresnel_f0(surf->ior, surf->albedo, surf->metallic);
//     float cosTheta = fmaxf(glms_vec3_dot(glms_vec3_scale(V, -1.0f), N), 0.0f);
//     vec3s base_F = schlick_fresnel(cosTheta, f0);

//     // Compute perfect mirror reflection direction.
//     vec3s base_direction = glms_vec3_normalize(glms_vec3_reflect(V, N));

//     vec3s origin = glms_vec3_add(q->hit_location, glms_vec3_scale(N, 0.001f));
//     ray base_ray = { origin, base_direction };

//     float gard_rotation = 0.1f * (float)sw_frame_index();

//     int32_t num_rays = min_reflected_rays + ((max_reflected_rays-min_reflected_rays) * surf->roughness);
//     ray perturbed_rays[max_reflected_rays];
//     gard_distribution_rays(base_ray, q->hit_surface->roughness, num_rays, gard_rotation, perturbed_rays);

//     for(int32_t i = 0; i < num_rays; i++)
//     {
//         // Perturb the reflection direction by mixing in a random direction from the hemisphere.
//         vec3s random_dir = random_in_hemisphere(N);
//         vec3s random_perturb_dir = glms_vec3_normalize(glms_vec3_mix(base_direction, random_dir, surf->roughness));
//         float gard_ratio = 0.75f;

//         perturbed_rays[i].direction = glms_vec3_normalize(glms_vec3_mix(random_perturb_dir, perturbed_rays[i].direction, gard_ratio));

//         out_rays[i].ray = perturbed_rays[i];
//         out_rays[i].weight = glms_vec3_scale(base_F, 1.0f / (float)num_rays);
//     }

//     return num_rays;
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

        //return glms_vec3_scale(albedo, cfg->camera->whitepoint);
    }
    else 
    {
        albedo = hit_surf->albedo;
    }

    // albedo = (vec3s){1.0f, 1.0f, 1.0f};

    float roughness;
    if(hit_surf->roughness_texture)
    {
        roughness = texture_sample_uv(hit_surf->roughness_texture, uv).r;
    }
    else 
    {
        roughness = hit_surf->roughness;
    }

    //roughness = powf(roughness, 1.3333f);
    //roughness = roughness * 0.5f;
    roughness = glm_clamp(roughness, 0.01, 1.0f);

    // roughness = 0.04f;

    float metallic;
    if(hit_surf->metallic_texture)
    {
        metallic = texture_sample_uv(hit_surf->metallic_texture, uv).r;
    }
    else 
    {
        metallic = hit_surf->metallic;
    }
    // metallic= 1.0f;

    vec3s N, T, B;
    T = glms_vec3_normalize(q->hit_tangent);
    B = glms_vec3_normalize(q->hit_bitangent);
    N = glms_vec3_normalize(q->hit_normal);

    if(hit_surf->normal_texture)
    {
        vec3s normal_map = glms_vec4_copy3(texture_sample_uv(hit_surf->normal_texture, uv));
        normal_map = glms_vec3_subs(glms_vec3_scale(normal_map, 2.0f), 1.0f);
        //normal_map.y = -normal_map.y;
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
        //mapped_normal.y = -mapped_normal.y;

        //float ndiff = 1.0f - fmaxf(0.0f, glms_vec3_dot(mapped_normal, N));
        //return (vec3s){ndiff, ndiff, ndiff};

        // Normal and view vectors.
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

    // Compute perfect mirror reflection direction.
    vec3s base_direction = glms_vec3_normalize(glms_vec3_reflect(V, N));

    vec3s origin = glms_vec3_add(q->hit_location, glms_vec3_scale(N, 0.001f));
    ray base_reflection_ray = { origin, base_direction };
    ray base_diffuse_ray = { origin, N };

    float gard_rotation = ((float)thread->sample_index / (float)cfg->num_samples) * (GLM_PI * 2.0);

    int32_t num_reflected_rays = cfg->min_reflected_rays + ((cfg->max_reflected_rays - cfg->min_reflected_rays) * roughness);
 
    // Roughness -> alpha
    float alpha = roughness * roughness;
    vec3s specular_F = glms_vec3_zero();
    for(int32_t i = 0; i < num_reflected_rays; i++)
    {
        // 1) Sample half-vector H from GGX
        vec3s H_local = sample_GGX_HalfVector(thread, alpha);
        vec3s H = to_world_space(H_local, N, T, B);  // Transform from local to world space

        // 2) Reflect the view vector about H
        vec3s L = glms_vec3_reflect(Vneg, H);

        float straightener = 1.0 - powf(roughness, 0.95f);
        L = glms_vec3_mix(L, base_direction, straightener);
        L = glms_vec3_normalize(L);

        // if(surf->roughness < 0.5)
        // {
        //     L = base_direction;
        // }

        float NdotL = fmaxf(glms_vec3_dot(N, L), 0.0f);
        if(NdotL <= 0.0f) {
            // This direction goes below the surface, discard or skip
            continue;
        }

        // 3) Compute the microfacet terms
        float NdotH = fmaxf(glms_vec3_dot(N, H), 0.0f);
        float VdotH = fmaxf(glms_vec3_dot(Vneg, H), 0.0f);

        // D: GGX normal distribution
        float denom = (NdotH * NdotH * (alpha - 1.0f) + 1.0f);
        denom = denom * denom; // squared
        float D = alpha / fmaxf(GLM_PI * denom, 0.0001f);

        // G: Geometry function (Smith with Schlick-GGX)
        float k = (roughness + 1.0f);
        k = (k * k) / 8.0f;
        float G1_V = NdotV / fmaxf((NdotV * (1.0f - k) + k), 0.0001f);
        float G1_L = NdotL / fmaxf((NdotL * (1.0f - k) + k), 0.0001f);
        float G = G1_V * G1_L;

        // F: Fresnel (Schlick)
        vec3s F = schlick_fresnel(VdotH, f0);
        specular_F = glms_vec3_add(specular_F, F);

        // 4) Full BRDF: f_spec(L,V) = (D * G * F) / (4*NdotV*NdotL)
        // We'll compute this now, but we also need the PDF
        // so we can do finalWeight = f_spec * (NdotL / PDF).
        // Because in path tracing: contribution = BRDF * cos(theta) / PDF
        // with cos(theta) = NdotL.

        // 5) PDF for GGX half-vector sampling
        //   pdf(L) = D(H)*NdotH / (4* VdotH)
        float pdf = (D * NdotH) / fmaxf(4.0f * VdotH, 0.0001f);

        // 6) Final sample contribution:
        //   f_spec(L,V) * (NdotL / pdf)
        //   but let's just do it in steps:
        float denomBRDF = 4.0f * NdotV * NdotL;
        vec3s specBRDF = glms_vec3_scale(F, (D * G) / fmaxf(denomBRDF, 0.0001f));

        // NdotL / pdf
        float ratio = (NdotL) / fmaxf(pdf, 0.0001f);
        vec3s reflection_weight = glms_vec3_scale(specBRDF, ratio);

        //reflection_weight = specBRDF;

        // 7) Push the new ray
        // The direction is L, the weight is reflection_weight / #samples
        // so we average over multiple samples.
        // (Or you can incorporate the "1 / num_specular_samples" into reflection_weight above.)
        vec3s final_weight = glms_vec3_scale(reflection_weight, 1.0f / (float)num_reflected_rays);

        // Create the new ray slightly offset along N
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




/*
 * Computes the refracted (transmitted) ray using Snell's law.
 * The function detects whether the ray is entering or exiting the material and calculates the
 * correct refraction direction. If total internal reflection occurs, the function falls back to
 * the reflection model.
 *
 * The Fresnel term (using Schlick's approximation) is used to determine the transmitted fraction
 * (i.e. 1 - F), which is returned as the weight.
//  */
// weighted_ray refraction_snell(query *q)
// {
//     surface_proxy *surf = q->hit_surface;
//     vec3s N = glms_vec3_normalize(q->hit_normal);
//     vec3s I = glms_vec3_normalize(q->ray.direction);

//     // Determine whether the ray is entering or exiting the material.
//     float cosi = fmaxf(-1.0f, fminf(1.0f, glms_vec3_dot(I, N)));
//     float etai = 1.0f;
//     float etat = surf->ior;

//     vec3s n = N;
//     if (cosi < 0.0f)
//     {
//         // Ray is entering.
//         cosi = -cosi;
//     }
//     else
//     {
//         // Ray is exiting; invert the normal and swap the indices.
//         n = glms_vec3_scale(N, -1.0f);
//         float temp = etai;
//         etai = etat;
//         etat = temp;
//     }
//     float eta = etai / etat;
//     float k = 1.0f - eta * eta * (1.0f - cosi * cosi);
//     vec3s refractDir;
//     if (k < 0.0f)
//     {
//         // Total internal reflection; fall back to microfacet reflection.
//         //return reflection_microfacet(q); @todo
//     }
//     else
//     {
//         refractDir = glms_vec3_add(glms_vec3_scale(I, eta),
//                                     glms_vec3_scale(n, (eta * cosi - sqrtf(k))));
//         refractDir = glms_vec3_normalize(refractDir);
//     }

//     // Compute Fresnel reflectance for a dielectric (F0 is about 0.04).
//     //vec3s F0 = { 0.04f, 0.04f, 0.04f };
//     vec3s f0 = fresnel_f0(surf->ior, surf->albedo, surf->metallic);
//     float cosTheta = fmaxf(glms_vec3_dot(glms_vec3_scale(I, -1.0f), N), 0.0f);
//     vec3s F = schlick_fresnel(cosTheta, f0);
//     // The transmitted fraction is (1 - F) per channel.
//     vec3s one = { 1.0f, 1.0f, 1.0f };
//     vec3s transmittance = glms_vec3_sub(one, F);

//     // Create the new refracted ray, offset slightly.
//     vec3s origin = glms_vec3_add(q->hit_location, glms_vec3_scale(refractDir, 0.001f));
//     ray newRay = { origin, refractDir };

//     weighted_ray ret;
//     ret.ray = newRay;
//     ret.weight = transmittance;
//     return ret;
// }


// vec3s environment_simple(ray ray)
// {
//     // Normalize the ray direction.
//     vec3s direction = glms_vec3_normalize(ray.direction);

//     // Map y-component of the ray direction to [0, 1].
//     float t = fmaxf(0.0f, fminf(1.0f, direction.y));

//     // Sky gradient: Interpolate between horizon color and zenith color.
//     vec3s horizonColor = glms_vec3_scale((vec3s) { 0.7f, 0.8f, 1.0f }, 2.0f); // Light blue.
//     vec3s zenithColor = { 0.0f, 0.05f, 0.2f }; // Deep blue.
//     vec3s skyColor = glms_vec3_mix(horizonColor, zenithColor, powf(t, 1.15f));

//     // Add a simple sun using a dot product to determine its intensity.
//     vec3s sunDirection   = glms_vec3_normalize((vec3s){ -1.0f, -1.0f, -1.0f });
//     float sunIntensity = powf(fmaxf(0.0f, glms_vec3_dot(direction, sunDirection)), 50.0f); // Sharpness of the sun.
//     vec3s sunColor = { 1.0f, 0.9f, 0.7f }; // Warm yellowish sun.
//     vec3s sunContribution = glms_vec3_scale(sunColor, sunIntensity* 10.0f);

//     // Combine the sky gradient and the sun contribution.
//     return glms_vec3_add(skyColor, sunContribution);
// }


// vec3s irradiance_simple(ray ray)
// {
//     // Normalize the ray direction.
//     vec3s direction = glms_vec3_normalize(ray.direction);

//     // Map y-component of the ray direction to [0, 1].
//     float t = fmaxf(0.0f, fminf(1.0f, direction.y));

//     t = 0.25f + t * 0.5f;

//     // Sky gradient: Interpolate between horizon color and zenith color.
//     vec3s horizonColor = glms_vec3_scale((vec3s) { 0.7f, 0.8f, 1.0f }, 2.0f); // Light blue.
//     vec3s zenithColor = { 0.0f, 0.05f, 0.2f }; // Deep blue.
//     vec3s skyColor = glms_vec3_mix(horizonColor, zenithColor, t);


//     // Combine the sky gradient and the sun contribution.
//     return  skyColor;
// }