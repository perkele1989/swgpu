
#pragma once 

#include <core.h>
#include <utils.h>

#include "common.h"

typedef struct vertex 
{
    vec3s position;
    vec3s normal;
    vec3s tangent;
    vec3s bitangent;
    vec2s uv;
} vertex;

typedef struct triangle 
{
    vertex v0, v1, v2;
    vec3s normal;
} triangle;

bool ray_triangle_intersection(ray r, triangle *t, float *out_distance);
vec3s triangle_barycentric(triangle *t, vec3s intersection_point);
vertex triangle_barycentric_interpolation(triangle *t, vec3s barycentric);
vec3s triangle_centroid(triangle *tri);
aabb aabb_from_triangles(triangle* tris, uint64_t num_tris);

typedef struct mesh_node
{
    aabb bounds;

    bool is_leaf;

    struct mesh_node* left;// only for non-leaf nodes
    struct mesh_node *right; // only for non-leaf nodes

    triangle *triangles; // only for leaf nodes
    uint32_t num_triangles; // only for leaf nodes
} mesh_node;

/** a mesh built up as a BVH */
typedef struct mesh 
{
    /** non-owning pointer to the root mesh */
    mesh_node *root;

    /** array of owning nodes */
    mesh_node* nodes;
    uint32_t num_nodes;

    /** array of owning triangles*/
    triangle* triangles;
    uint32_t num_triangles;
} mesh;


/** loads the mesh from the given path, will build a BVH for it  */
bool mesh_load_fbx(mesh *out_mesh, char const *path);

void mesh_destroy(mesh* m);

bool mesh_traverse(mesh_node *node, ray ray, float *closest_hit_distance, triangle **hit_triangle);

mesh_node* mesh_build_bvh_node(mesh* mesh, triangle* triangles, uint64_t num_triangles, uint64_t treshold);