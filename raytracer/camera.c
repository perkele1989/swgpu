#include "camera.h"
#include "common.h"

void camera_apply(camera *cam)
{
    cam->derivatives.projection_matrix = glms_perspective(cam->vertical_fov, cam->aspect, cam->near, cam->far);
    cam->derivatives.inv_projection_matrix = glms_mat4_inv(cam->derivatives.projection_matrix);

    mat4s rotation_matrix = glms_quat_mat4(cam->rotation);
    mat4s translation_matrix = glms_translate_make(cam->position);

    cam->derivatives.transform_matrix = glms_mat4_mul( translation_matrix, rotation_matrix);
    cam->derivatives.view_matrix = glms_mat4_inv(cam->derivatives.transform_matrix); 

    cam->derivatives.forward = glms_mat4_mulv3(cam->derivatives.transform_matrix, world_forward, 0.0f);
    cam->derivatives.right = glms_mat4_mulv3(cam->derivatives.transform_matrix, world_right, 0.0f);
    cam->derivatives.up = glms_mat4_mulv3(cam->derivatives.transform_matrix, world_up, 0.0f);
}


void camera_init(camera *cam, float vfov, float aspect, float near, float far)
{
    cam->position = (vec3s){0.0f, 0.0f, 0.0f};
    cam->rotation = glms_quat_identity();

    cam->vertical_fov = vfov;
    cam->aspect = aspect;
    cam->near = near;
    cam->far = far;

    camera_apply(cam);
}

vec3s camera_view_vector(camera *cam, float ndc_x, float ndc_y)
{
    vec4s dir_camera = {ndc_x, ndc_y, -1.0f, 1.0f};
    vec4s dir_camera_space = glms_mat4_mulv(cam->derivatives.inv_projection_matrix, dir_camera);
    if(dir_camera_space.w != 0.0f)
    {
        dir_camera_space.x /= dir_camera_space.w;
        dir_camera_space.y /= dir_camera_space.w;
        dir_camera_space.z /= dir_camera_space.w;
        dir_camera_space.w = 0.0f;
    }

    vec4s dir_world = glms_mat4_mulv(cam->derivatives.transform_matrix, dir_camera_space);
    vec3s dw = {dir_world.x, dir_world.y, dir_world.z};
    dw = glms_vec3_normalize(dw);

    return dw;
}
