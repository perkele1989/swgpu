
#include "scene.h"

#include "mesh.h"


void entity_set_transform(entity *m, vec3s position, versors rotation)
{
    mat4s mat_t = glms_translate_make(position);
    mat4s mat_r = glms_quat_mat4(rotation);

    m->transform = glms_mat4_mul(mat_t, mat_r);
    m->inverse_transform = glms_mat4_inv(m->transform);
    m->normal = glms_mat3_transpose(glms_mat4_pick3(m->inverse_transform));
    m->inverse_normal = glms_mat3_inv(m->normal);
}

void scene_create(scene *s)
{
    array_create(&s->entities, sizeof(entity), 8);
}

query scene_query(scene *scene, ray in_ray)
{
    query returner;
    returner.ray = in_ray;
    returner.scene = scene;
    returner.hit = false;
    returner.hit_distance = FLT_MAX;
    returner.hit_location = glms_vec3_zero();
    returner.hit_normal = world_up;
    returner.hit_tangent = world_right;
    returner.hit_bitangent = world_forward;
    returner.hit_uv = (vec2s){0.0f, 0.0f};
    returner.hit_surface = NULL;

    entity *hit_entity = NULL;

    for (int32_t m = 0; m < scene->entities.size; m++)
    {
        entity *current_entity = array_at(&scene->entities, m);
        float mesh_hit_distance = FLT_MAX;
        triangle *hit_triangle = NULL;

        ray mesh_ray = in_ray;
        mesh_ray.origin = glms_mat4_mulv3(current_entity->inverse_transform, mesh_ray.origin, 1.0f);
        mesh_ray.direction = glms_vec3_normalize(glms_mat3_mulv(current_entity->inverse_normal, mesh_ray.direction));

        // traverse mesh bvh
        if (mesh_traverse(current_entity->mesh->root, mesh_ray, &mesh_hit_distance, &hit_triangle))
        {
            if (mesh_hit_distance < returner.hit_distance)
            {
                hit_entity = current_entity;
                returner.hit = true;
                returner.hit_distance = mesh_hit_distance;
                returner.hit_location = glms_vec3_add(mesh_ray.origin, glms_vec3_scale(mesh_ray.direction, mesh_hit_distance));
                vec3s bary = triangle_barycentric(hit_triangle, returner.hit_location);
                vertex interp_vert = triangle_barycentric_interpolation(hit_triangle, bary);
                returner.hit_normal = interp_vert.normal;
                returner.hit_tangent = interp_vert.tangent;
                returner.hit_bitangent = interp_vert.bitangent;
                returner.hit_uv = interp_vert.uv;
                returner.hit_surface = hit_entity->surface;
            }
        }
    }

    if(hit_entity)
    {
        returner.hit_location = glms_mat4_mulv3(hit_entity->transform, returner.hit_location, 1.0f);
        returner.hit_normal = glms_vec3_normalize(glms_mat3_mulv(hit_entity->normal, returner.hit_normal));
        returner.hit_tangent =  glms_vec3_normalize(glms_mat3_mulv(hit_entity->normal, returner.hit_tangent));
        returner.hit_bitangent =  glms_vec3_normalize(glms_mat3_mulv(hit_entity->normal, returner.hit_bitangent));
    }

    return returner;
}






fast_query scene_fast_query(scene *scene, ray in_ray)
{
    fast_query returner;
    returner.ray = in_ray;
    returner.scene = scene;
    returner.hit = false;

    float hit_distance = FLT_MAX;

    entity *hit_entity = NULL;

    for (int32_t m = 0; m < scene->entities.size; m++)
    {
        entity *current_entity = array_at(&scene->entities, m);
        float mesh_hit_distance = FLT_MAX;
        triangle *hit_triangle = NULL;

        ray mesh_ray = in_ray;
        mesh_ray.origin = glms_mat4_mulv3(current_entity->inverse_transform, mesh_ray.origin, 1.0f);
        mesh_ray.direction = glms_vec3_normalize(glms_mat3_mulv(current_entity->inverse_normal, mesh_ray.direction));

        // traverse mesh bvh
        if (mesh_traverse(current_entity->mesh->root, mesh_ray, &mesh_hit_distance, &hit_triangle))
        {
            if (mesh_hit_distance < hit_distance)
            {
                hit_entity = current_entity;
                returner.hit = true;
                hit_distance = mesh_hit_distance;
            }
        }
    }

    return returner;
}



void scene_destroy(scene *s)
{
    array_destroy(&s->entities);
}
