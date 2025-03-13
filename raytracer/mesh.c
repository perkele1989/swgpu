
#include "mesh.h"

#include "ufbx.h"


vec2s ufbx2_to_glms(ufbx_vec2 in)
{
    return (vec2s){in.x, in.y};
}

vec3s ufbx3_to_glms(ufbx_vec3 in)
{
    return (vec3s){in.x, in.y, in.z};
}

vec4s ufbx4_to_glms(ufbx_vec4 in)
{
    return (vec4s){in.x, in.y, in.z, in.w};
}



bool mesh_load_fbx(mesh *out_mesh, char const *path)
{
    ufbx_load_opts opts = { 0 }; // Optional, pass NULL for defaults
    opts.obj_merge_groups = true;
    opts.obj_merge_objects = true;
    opts.generate_missing_normals = true;
    opts.normalize_normals = true;
    opts.normalize_tangents = true;
    opts.geometry_transform_handling = UFBX_GEOMETRY_TRANSFORM_HANDLING_MODIFY_GEOMETRY;

    ufbx_error error; // Optional, pass NULL if you don't care about errors
    ufbx_scene *scene = ufbx_load_file(path, &opts, &error);
    if (!scene) 
    {
        printf("Failed to load: %s\n", error.description.data);
        return false;
    }

    if(scene->meshes.count == 0)
    {
        printf("zero meshes in file\n");
        return false;
    }

    if(scene->meshes.count > 1)
    {
        printf("more than one mesh in fbx, we will simply use the first one\n");
    }

    ufbx_mesh* in_mesh = scene->meshes.data[0];

    if(!in_mesh->vertex_position.exists)
    {
        printf("no vertex positions for mesh!\n");
        return false;
    }

    if(!in_mesh->vertex_normal.values.data)
    {
        printf("no vertex normals for mesh!\n");
        return false;
    }

    if(!in_mesh->vertex_tangent.values.data)
    {
        printf("no vertex tangents for mesh!\n");
        return false;
    }

    if(!in_mesh->vertex_bitangent.values.data)
    {
        printf("no vertex bitangents for mesh!\n");
        return false;
    }

    bool warn_faces = false;

    out_mesh->triangles = malloc(sizeof(triangle) * in_mesh->faces.count);
    out_mesh->num_triangles = 0;
    memset(out_mesh->triangles, 0, sizeof(triangle) * in_mesh->faces.count);

    out_mesh->nodes =  malloc(sizeof(mesh_node) * in_mesh->faces.count * 2);
    out_mesh->num_nodes = 0;
    memset(out_mesh->nodes, 0, sizeof(mesh_node) * in_mesh->faces.count * 2);


    for(size_t f = 0; f < in_mesh->faces.count; f++)
    {
        ufbx_face* face = &in_mesh->faces.data[f];
        if(face->num_indices != 3)
        {
            if(!warn_faces)
            {
                printf("model is not triangulated, will ignore faces.\n");
                warn_faces = true;
            }
            continue;
        }

        vec3s pos0 = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_position, face->index_begin + 0));
        vec3s pos1 = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_position, face->index_begin + 1));
        vec3s pos2 = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_position, face->index_begin + 2));

        triangle new_triangle;
        new_triangle.v0.position = pos0;
        new_triangle.v1.position = pos1;
        new_triangle.v2.position = pos2;

        new_triangle.v0.normal = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_normal, face->index_begin + 0));
        new_triangle.v1.normal = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_normal, face->index_begin + 1));
        new_triangle.v2.normal = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_normal, face->index_begin + 2));

        new_triangle.v0.tangent = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_tangent, face->index_begin + 0));
        new_triangle.v1.tangent = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_tangent, face->index_begin + 1));
        new_triangle.v2.tangent = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_tangent, face->index_begin + 2));

        new_triangle.v0.bitangent = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_bitangent, face->index_begin + 0));
        new_triangle.v1.bitangent = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_bitangent, face->index_begin + 1));
        new_triangle.v2.bitangent = ufbx3_to_glms(ufbx_get_vertex_vec3(&in_mesh->vertex_bitangent, face->index_begin + 2));

        if(in_mesh->vertex_uv.values.data)
        {
            new_triangle.v0.uv = ufbx2_to_glms(ufbx_get_vertex_vec2(&in_mesh->vertex_uv, face->index_begin + 0));
            new_triangle.v1.uv = ufbx2_to_glms(ufbx_get_vertex_vec2(&in_mesh->vertex_uv, face->index_begin + 1));
            new_triangle.v2.uv = ufbx2_to_glms(ufbx_get_vertex_vec2(&in_mesh->vertex_uv, face->index_begin + 2));
        }
        else 
        {
            new_triangle.v0.uv = glms_vec2_zero();
            new_triangle.v1.uv = glms_vec2_zero();
            new_triangle.v2.uv = glms_vec2_zero();
        }

        vec3s edge1 = glms_vec3_normalize(glms_vec3_sub(new_triangle.v1.position, new_triangle.v0.position));
        vec3s edge2 = glms_vec3_normalize(glms_vec3_sub(new_triangle.v2.position, new_triangle.v0.position));
        new_triangle.normal = glms_vec3_normalize(glms_vec3_cross(edge1, edge2));
        out_mesh->triangles[out_mesh->num_triangles++] = new_triangle;
    }

    ufbx_free_scene(scene);

    out_mesh->root = mesh_build_bvh_node(out_mesh, out_mesh->triangles, out_mesh->num_triangles, 8);
    if(!out_mesh->root)
    {
        mesh_destroy(out_mesh);
        return false;
    }

    return true;
}


void mesh_destroy(mesh* out_mesh)
{
    free(out_mesh->nodes);
    free(out_mesh->triangles);
}


int32_t _triangle_axis_cmp(axis* axis, triangle *a, triangle* b)
{
    vec3s centroid_a = triangle_centroid(a);
    vec3s centroid_b = triangle_centroid(b);
    if(*axis == axis_x)
    {
        return centroid_a.x == centroid_b.x ? 0 : centroid_a.x < centroid_b.x ? -1 : 1;
    }
    if(*axis == axis_y)
    {
        return centroid_a.y == centroid_b.y ? 0 : centroid_a.y < centroid_b.y ? -1 : 1;
    }

    return centroid_a.z == centroid_b.z ? 0 : centroid_a.z < centroid_b.z ? -1 : 1;
}



mesh_node* mesh_build_bvh_node(mesh* mesh, triangle* triangles, uint64_t num_triangles, uint64_t treshold)
{
    if(num_triangles <= 1)
    {
        printf("1 or less triangles, error! treshold must be at least 2. returning null");
        return NULL;
    }

    mesh_node *new_node = mesh->nodes + (mesh->num_nodes++);
    new_node->bounds = aabb_from_triangles(triangles, num_triangles);

    if(num_triangles <= treshold)
    {
        new_node->is_leaf = true;
        new_node->num_triangles = num_triangles;
        new_node->triangles = triangles;

        return new_node;
    }

    new_node->is_leaf = false;

    axis longest_axis = aabb_longest_axis(new_node->bounds);

    qsort_s(triangles, num_triangles, sizeof(triangle), _triangle_axis_cmp, &longest_axis);

    uint64_t mid = num_triangles / 2;
    
    new_node->left = mesh_build_bvh_node(mesh, triangles, mid, treshold);
    new_node->right = mesh_build_bvh_node(mesh, triangles + mid, num_triangles - mid, treshold);

    return new_node;
}



bool mesh_traverse(mesh_node *node, ray ray, float *closest_hit_distance, triangle **hit_triangle)
{
    if (!node)
        return false;

    // Check ray-AABB intersection
    float aabb_distance;
    if (!ray_aabb_intersection(node->bounds,ray, &aabb_distance))
        return false;

    bool hit = false;

    if (node->is_leaf)
    {
        // Test all triangles in the leaf node
        for (size_t i = 0; i < node->num_triangles; i++)
        {
            triangle *tri = &node->triangles[i];
            float tri_distance;
            if (ray_triangle_intersection(ray, tri, &tri_distance) && tri_distance < *closest_hit_distance)
            {
                *closest_hit_distance = tri_distance;
                *hit_triangle = tri;
                hit = true;
            }
        }
    }
    else
    {
        // Recurse into child nodes
        bool left_hit = mesh_traverse(node->left, ray, closest_hit_distance, hit_triangle);
        bool right_hit = mesh_traverse(node->right, ray, closest_hit_distance, hit_triangle);
        hit = left_hit || right_hit;
    }

    return hit;
}


bool ray_triangle_intersection(ray r, triangle *t, float *out_distance)
{
    *out_distance = 0.0f;

    //if(glms_vec3_dot(r.direction, t->normal) > 0.0f)
    //   return false;

    // Check for intersection with the plane of the triangle
    plane p = { .point = t->v0.position, .normal = t->normal };

    float plane_distance = 0.0f;
    if (!ray_plane_intersection(r, p, &plane_distance))
    {
        return false;
    }

    // Compute the intersection point
    vec3s intersection_point = glms_vec3_add(r.origin, glms_vec3_scale(r.direction, plane_distance));

    // Perform inside-outside test for the triangle

    // Edge 0
    vec3s edge0 = glms_vec3_sub(t->v1.position, t->v0.position);
    vec3s vp0 = glms_vec3_sub(intersection_point, t->v0.position);
    vec3s c0 = glms_vec3_cross(edge0, vp0);
    if (glms_vec3_dot(t->normal, c0) < 0.0f)
        return false;

    // Edge 1
    vec3s edge1_check = glms_vec3_sub(t->v2.position, t->v1.position);
    vec3s vp1 = glms_vec3_sub(intersection_point, t->v1.position);
    vec3s c1 = glms_vec3_cross(edge1_check, vp1);
    if (glms_vec3_dot(t->normal, c1) < 0.0f)
        return false;

    // Edge 2
    vec3s edge2_check = glms_vec3_sub(t->v0.position, t->v2.position);
    vec3s vp2 = glms_vec3_sub(intersection_point, t->v2.position);
    vec3s c2 = glms_vec3_cross(edge2_check, vp2);
    if (glms_vec3_dot(t->normal, c2) < 0.0f)
        return false;

    // If we passed all edge checks, the ray intersects the triangle
    *out_distance = plane_distance;

    return true;
}



vec3s triangle_barycentric(triangle *t, vec3s p)
{
    // Triangle vertices
    vec3s a = t->v0.position;
    vec3s b = t->v1.position;
    vec3s c = t->v2.position;

    vec3s v0 = glms_vec3_sub(b, a);
    vec3s v1 = glms_vec3_sub(c, a);
    vec3s v2 = glms_vec3_sub(p, a);

    float d00 = glms_vec3_dot(v0, v0);
    float d01 = glms_vec3_dot(v0, v1);
    float d11 = glms_vec3_dot(v1, v1);
    float d20 = glms_vec3_dot(v2, v0);
    float d21 = glms_vec3_dot(v2, v1);

    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    return (vec3s){u, v, w};

}

vec3s triangle_barycentric2(triangle *t, vec3s intersection_point)
{
    // Triangle vertices
    vec3s v0 = t->v0.position;
    vec3s v1 = t->v1.position;
    vec3s v2 = t->v2.position;

    // Compute the full triangle's normal
    vec3s edge1 = glms_vec3_sub(v1, v0);
    vec3s edge2 = glms_vec3_sub(v2, v0);
    vec3s triangle_normal = glms_vec3_cross(edge1, edge2);
    float triangle_area_inv = 1.0f / glms_vec3_norm(triangle_normal);

    // Sub-triangle areas
    vec3s p_to_v0 = glms_vec3_sub(intersection_point, v0);
    vec3s p_to_v1 = glms_vec3_sub(intersection_point, v1);
    vec3s p_to_v2 = glms_vec3_sub(intersection_point, v2);

    vec3s normal_u = glms_vec3_cross(glms_vec3_sub(v2, v1), p_to_v1); // (v1, v2, intersection_point)
    vec3s normal_v = glms_vec3_cross(glms_vec3_sub(v0, v2), p_to_v2); // (v2, v0, intersection_point)
    vec3s normal_w = glms_vec3_cross(glms_vec3_sub(v1, v0), p_to_v0); // (v0, v1, intersection_point)

    float u = glms_vec3_dot(normal_u, triangle_normal) * triangle_area_inv;
    float v = glms_vec3_dot(normal_v, triangle_normal) * triangle_area_inv;
    float w = glms_vec3_dot(normal_w, triangle_normal) * triangle_area_inv;

    return (vec3s){{u, v, w}};
}

vertex triangle_barycentric_interpolation(triangle *t, vec3s barycentric)
{
    vertex returner;
    returner.position = glms_vec3_zero();
    returner.position = glms_vec3_add(returner.position, glms_vec3_scale(t->v0.position, barycentric.x));
    returner.position = glms_vec3_add(returner.position, glms_vec3_scale(t->v1.position, barycentric.y));
    returner.position = glms_vec3_add(returner.position, glms_vec3_scale(t->v2.position, barycentric.z));

    returner.normal = glms_vec3_zero();
    returner.normal = glms_vec3_add(returner.normal, glms_vec3_scale(t->v0.normal, barycentric.x));
    returner.normal = glms_vec3_add(returner.normal, glms_vec3_scale(t->v1.normal, barycentric.y));
    returner.normal = glms_vec3_add(returner.normal, glms_vec3_scale(t->v2.normal, barycentric.z));
    returner.normal = glms_vec3_normalize(returner.normal);

    returner.tangent = glms_vec3_zero();
    returner.tangent = glms_vec3_add(returner.tangent, glms_vec3_scale(t->v0.tangent, barycentric.x));
    returner.tangent = glms_vec3_add(returner.tangent, glms_vec3_scale(t->v1.tangent, barycentric.y));
    returner.tangent = glms_vec3_add(returner.tangent, glms_vec3_scale(t->v2.tangent, barycentric.z));
    returner.tangent = glms_vec3_normalize(returner.tangent);

    returner.bitangent = glms_vec3_zero();
    returner.bitangent = glms_vec3_add(returner.bitangent, glms_vec3_scale(t->v0.bitangent, barycentric.x));
    returner.bitangent = glms_vec3_add(returner.bitangent, glms_vec3_scale(t->v1.bitangent, barycentric.y));
    returner.bitangent = glms_vec3_add(returner.bitangent, glms_vec3_scale(t->v2.bitangent, barycentric.z));
    returner.bitangent = glms_vec3_normalize(returner.bitangent);

    returner.uv = glms_vec2_zero();
    returner.uv = glms_vec2_add(returner.uv, glms_vec2_scale(t->v0.uv, barycentric.x));
    returner.uv = glms_vec2_add(returner.uv, glms_vec2_scale(t->v1.uv, barycentric.y));
    returner.uv = glms_vec2_add(returner.uv, glms_vec2_scale(t->v2.uv, barycentric.z));

    return returner;
}



aabb aabb_from_triangles(triangle* tris, uint64_t num_tris)
{
    aabb returner = aabb_empty();

    for(uint64_t i = 0; i < num_tris; i++)
    {
        aabb_add_point(&returner, tris->v0.position);
        aabb_add_point(&returner, tris->v1.position);
        aabb_add_point(&returner, tris->v2.position);
        tris++;
    }
        

    return returner;
}

vec3s triangle_centroid(triangle *tri)
{
    return glms_vec3_scale(glms_vec3_add(glms_vec3_add(tri->v0.position, tri->v1.position), tri->v2.position), 1.0f / 3.0f);
}