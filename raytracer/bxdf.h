
#pragma once 

#include "pathtracer.h"

vec3s fresnel_f0(float ior, vec3s albedo, float metallic);

// Schlick's approximation for Fresnel reflectance.
// Returns: F = F0 + (1 - F0) * (1 - cosTheta)^5
vec3s schlick_fresnel(float cosTheta, vec3s f0);

// Returns a random unit vector in the hemisphere around the given normal.
// (A very basic implementation; in production code you’d want a better sampler.)
vec3s random_in_hemisphere(tracer_thread* thread, vec3s normal);

vec3s material_light(tracer_thread* thread, query* q, vec3s weight, int32_t depth);

vec3s material_default(tracer_thread* thread, query* q, vec3s weight, int32_t depth);

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
// vec3s light_cooktorrance(query *q);

/*
 * Computes the reflected ray using a microfacet approach.
 * The perfect mirror reflection is computed first, and then it is perturbed by a random
 * offset that is proportional to the surface roughness. (A rough surface produces a blurrier
 * reflection.) The Fresnel term (using Schlick's approximation) is used as the reflectance weight.
 */
// int32_t reflection_microfacet(query *q, weighted_ray* out_rays);
// int32_t reflection_microfacet_multi_rand(query *q, weighted_ray* out_rays);
// int32_t reflection_microfacet_multi_gard(query *q, weighted_ray* out_rays);

/*
 * Computes the refracted (transmitted) ray using Snell's law.
 * The function detects whether the ray is entering or exiting the material and calculates the
 * correct refraction direction. If total internal reflection occurs, the function falls back to
 * the reflection model.
 *
 * The Fresnel term (using Schlick's approximation) is used to determine the transmitted fraction
 * (i.e. 1 - F), which is returned as the weight.
//  */
// weighted_ray refraction_snell(query *q);

// vec3s environment_simple(ray ray);
// vec3s irradiance_simple(ray ray);