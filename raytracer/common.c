
#include "common.h"


/** Golden Angle Ray Distribution */
void gard_distribution_rays(ray base_ray, float roughness, int num_rays, float frame_rotation, ray *out_rays)
{
    vec3s base_direction = base_ray.direction;

    // Create an orthonormal basis around the base_direction
    // z-polar method (seems best for now):
    vec3s tangent = (fabsf(base_direction.z) < 0.999999f) ? 
                    glms_vec3_cross(base_direction, (vec3s){{0.0f, 0.0f, 1.0f}}) : 
                    glms_vec3_cross(base_direction, (vec3s){{0.0f, 1.0f, 0.0f}});

    tangent = glms_vec3_normalize(tangent);
    vec3s bitangent = glms_vec3_cross(base_direction, tangent);

    // Scale the cone angle by roughness
    float cone_angle = roughness * M_PI_2; // Roughness scales to [0, 90 degrees]

    // Precompute useful constants
    float golden_angle = 2.39996322972865332f; // ~137.5 degrees in radians
    float cone_cos = cosf(cone_angle);

    for (int i = 0; i < num_rays; i++) {
        // Calculate seed positions using the sunflower pattern
        float radius = sqrtf((float)i / num_rays);    // Radius within the cone
        float theta = i * golden_angle + frame_rotation; // Apply frame rotation

        // Map radius to spherical coordinates
        float phi = radius * sinf(cone_angle); // Angle from the base direction
        float z = cone_cos + (1.0f - cone_cos) * radius; // Interpolate in the cone

        // Compute local offsets
        float x = sqrtf(1.0f - z * z) * cosf(theta);
        float y = sqrtf(1.0f - z * z) * sinf(theta);

        // Transform into world space using the tangent/bitangent frame
        vec3s direction = glms_vec3_add(
            glms_vec3_add(glms_vec3_scale(tangent, x), glms_vec3_scale(bitangent, y)),
            glms_vec3_scale(base_direction, z)
        );

        // Normalize to ensure unit vector
        out_rays[i].origin = base_ray.origin;
        out_rays[i].direction = glms_vec3_normalize(direction);
    }
}

bool ray_sphere_intersection(ray r, vec3s center, float radius, float *out_distance)
{
    *out_distance = 0.0f;
    // Vector from the ray origin to the sphere center
    vec3s L = glms_vec3_sub(center, r.origin);

    // Project L onto the ray direction
    float tca = glms_vec3_dot(L, r.direction);

    // Compute the squared distance from the sphere center to the ray
    float d2 = glms_vec3_dot(L, L) - tca * tca;

    // If the distance is greater than the sphere's radius, there's no hit
    float radius2 = radius * radius;
    if (d2 > radius2)
    {
        return false;
    }

    // Compute the distance from the closest point to the intersection point
    float thc = sqrtf(radius2 - d2);

    // Compute the two intersection points
    float t0 = tca - thc;
    float t1 = tca + thc;

    // If both intersections are behind the ray, there's no hit
    if (t0 < 0 && t1 < 0)
    {
        return false;
    }

    // Use the nearest positive intersection point
    float t = (t0 > 0) ? t0 : t1;
    *out_distance = t;

    return true;
}


aabb aabb_empty()
{
    aabb box;
    box.min.x = FLT_MAX;
    box.min.y = FLT_MAX;
    box.min.z = FLT_MAX;
    box.max.x = -FLT_MAX;
    box.max.y = -FLT_MAX;
    box.max.z = -FLT_MAX;
    return box;
}

void aabb_add_point(aabb *box, vec3s point)
{
    box->min.x = fminf(box->min.x, point.x);
    box->min.y = fminf(box->min.y, point.y);
    box->min.z = fminf(box->min.z, point.z);

    box->max.x = fmaxf(box->max.x, point.x);
    box->max.y = fmaxf(box->max.y, point.y);
    box->max.z = fmaxf(box->max.z, point.z);
}

// Optionally, merge two AABBs into one that encloses both.
aabb aabb_merge(aabb a, aabb b)
{
    aabb dst;
    dst.min.x = fminf(a.min.x, b.min.x);
    dst.min.y = fminf(a.min.y, b.min.y);
    dst.min.z = fminf(a.min.z, b.min.z);
    dst.max.x = fmaxf(a.max.x, b.max.x);
    dst.max.y = fmaxf(a.max.y, b.max.y);
    dst.max.z = fmaxf(a.max.z, b.max.z);
    return dst;
}

bool ray_plane_intersection(ray r, plane p, float *out_distance)
{
    *out_distance = 0.0f;

    // Dot product of ray direction and plane normal
    float denom = glms_vec3_dot(p.normal, r.direction);
    if (fabsf(denom) < 1e-6)
    {
        // The ray is parallel to the plane
        return false;
    }

    // Compute the distance along the ray to the intersection point
    vec3s diff = glms_vec3_sub(p.point, r.origin);
    float t = glms_vec3_dot(diff, p.normal) / denom;

    if (t < 0)
    { 
        // The intersection is behind the ray origin
        return false;
    }

    *out_distance = t;

    return true;
}

bool ray_aabb_intersection(aabb box, ray r, float *out_distance)
{
    vec3s inv_dir = glms_vec3_div((vec3s){{1.0f, 1.0f, 1.0f}}, r.direction); // Compute inverse direction

    float tmin, tmax, tymin, tymax, tzmin, tzmax;

    // Compute t-values for the x planes
    tmin = (box.min.x - r.origin.x) * inv_dir.x;
    tmax = (box.max.x - r.origin.x) * inv_dir.x;
    if (tmin > tmax) glm_swapf(&tmin, &tmax);

    // Compute t-values for the y planes
    tymin = (box.min.y - r.origin.y) * inv_dir.y;
    tymax = (box.max.y - r.origin.y) * inv_dir.y;
    if (tymin > tymax) glm_swapf(&tymin, &tymax);

    // Check for overlap in t intervals
    if ((tmin > tymax) || (tymin > tmax))
        return false;

    // Update tmin and tmax to include y intervals
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    // Compute t-values for the z planes
    tzmin = (box.min.z - r.origin.z) * inv_dir.z;
    tzmax = (box.max.z - r.origin.z) * inv_dir.z;
    if (tzmin > tzmax) glm_swapf(&tzmin, &tzmax);

    // Check for overlap in t intervals
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    // Update tmin and tmax to include z intervals
    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;

    // If `tmax` is positive, there's an intersection
    if (tmax >= 0.0f)
    {
        *out_distance = (tmin >= 0.0f) ? tmin : 0.0f; // Use `tmin` if positive; otherwise, use 0.0f for inside box
        return true;
    }

    return false; // No valid intersection
}

axis aabb_longest_axis(aabb box)
{
    vec3s box_size = glms_vec3_sub(box.max, box.min);
    if(box_size.x >= box_size.y && box_size.x >= box_size.z)
    {
        return axis_x;
    }

    if(box_size.y >= box_size.z)
    {
        return axis_y;
    }

    return axis_z;
}


// Convert an sRGB color (in [0,1]) to linear space.
vec3s srgb_to_linear(vec3s c) {
    vec3s result;
    result.x = (c.x <= 0.04045f) ? (c.x / 12.92f)
                                 : powf((c.x + 0.055f) / 1.055f, 2.4f);
    result.y = (c.y <= 0.04045f) ? (c.y / 12.92f)
                                 : powf((c.y + 0.055f) / 1.055f, 2.4f);
    result.z = (c.z <= 0.04045f) ? (c.z / 12.92f)
                                 : powf((c.z + 0.055f) / 1.055f, 2.4f);
    return result;
}

// Convert a linear color (in [0,1]) to sRGB space.
vec3s linear_to_srgb(vec3s c) {
    vec3s result;
    result.x = (c.x <= 0.0031308f) ? (c.x * 12.92f)
                                  : (1.055f * powf(c.x, 1.0f / 2.4f) - 0.055f);
    result.y = (c.y <= 0.0031308f) ? (c.y * 12.92f)
                                  : (1.055f * powf(c.y, 1.0f / 2.4f) - 0.055f);
    result.z = (c.z <= 0.0031308f) ? (c.z * 12.92f)
                                  : (1.055f * powf(c.z, 1.0f / 2.4f) - 0.055f);
    return result;
}