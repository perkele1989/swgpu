
#pragma once 

#include "core.h"

typedef struct camera 
{
    vec3s position;
    versors rotation;

    float vertical_fov;

    float aspect;
    float near;
    float far;

    float whitepoint;

    // derivatives of the above values, calculated on camera_apply()
    struct 
    {
        vec3s forward;
        vec3s right;
        vec3s up;

        mat4s view_matrix;
        mat4s projection_matrix;

        mat4s transform_matrix;
        mat4s inv_projection_matrix;
    } derivatives;
} camera;


void camera_apply(camera *cam);
void camera_init(camera *cam, float vfov, float aspect, float near, float far);
vec3s camera_view_vector(camera *cam, float ndc_x, float ndc_y);
