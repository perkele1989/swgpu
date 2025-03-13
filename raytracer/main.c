
#include "core.h"
#include "bxdf.h"
#include "camera.h"
#include "texture.h"
#include "scene.h"
#include "mesh.h"
#include "pathtracer.h"

#include <threads.h>




typedef struct app_state
{
    ivec2s render_size;

    float pitch;
    float yaw;
    camera camera;
    scene scene;

    surface surf_light;
    surface surf_gold[5]; 
    surface surf_plastic[5];

    texture cube_env;

    mesh mesh_bust;
    mesh mesh_cube;

    mesh mesh_beretta;
    texture tex_beretta_albedo;
    texture tex_beretta_roughness;
    texture tex_beretta_normal;

    surface surf_beretta;

    bool movement_enabled;

    tracer_cfg tracer_cfg;
    tracer tracer;

    uint32_t fps_frames;
    double fps_seconds;

    double seconds_since_preview;
} app_state; 


static app_state app;

void render_viewport(color *buffer)
{
    int32_t w = app.render_size.x/8;
    int32_t h = app.render_size.y/8;
    for(int32_t y = 0; y < h; y++)
    {
        for(int32_t x = 0; x < w; x++)
        {
            float u = (float)x / (float)w;
            float v = (float)y / (float)h;
        
            float ndc_x = u * 2.0f - 1.0f;
            float ndc_y = v * 2.0f - 1.0f;
        
            vec3s view_vector = camera_view_vector(&app.camera, ndc_x, ndc_y);
        
            ray initial_ray;
            initial_ray.direction = view_vector;
            initial_ray.origin = app.camera.position;
            
            vec3s output;
            fast_query q = scene_fast_query(&app.scene, initial_ray);
            if(q.hit)
            {
                output = glms_vec3_scale((vec3s){1.0f, 1.0f, 1.0f}, app.camera.whitepoint);
            }
            else 
            {
                //output = glms_vec3_scale(glms_vec4_copy3(texture_sample_uvw(app.scene.environment_hdr, initial_ray.direction)), app.scene.environment_intensity);
                output = (vec3s){0.2f, 0.2f, 0.2f};
            }

            output = tonemap_reinhard2(output, app.camera.whitepoint);
            buffer[x + y * app.render_size.x] = vec3_to_color(output);
        }
    }
}

void setup_assets()
{
    //texture_load(&app.cube_env, "../assets/symmetrical_garden_02_4k.hdr");
    //texture_load(&app.cube_env, "../assets/rostock_laage_airport_4k.hdr");
    texture_load(&app.cube_env, "../assets/the_sky_is_on_fire_4k.hdr");

    mesh_load_fbx(&app.mesh_bust, "../assets/bust.fbx");
    mesh_load_fbx(&app.mesh_cube, "../assets/cube2.fbx");
    mesh_load_fbx(&app.mesh_beretta, "../assets/beretta.fbx");

    texture_load(&app.tex_beretta_albedo, "../assets/beretta_albedo.png");
    texture_load(&app.tex_beretta_normal, "../assets/beretta_normal.png");
    texture_load(&app.tex_beretta_roughness, "../assets/beretta_roughness.png");

    app.surf_beretta.metallic = 0.0f;
    app.surf_beretta.emissive = glms_vec3_zero();
    app.surf_beretta.ior = 1.450f;
    app.surf_beretta.material_func = material_default;
    app.surf_beretta.roughness = 1.0f;
    app.surf_beretta.albedo_texture = &app.tex_beretta_albedo;
    app.surf_beretta.roughness_texture = &app.tex_beretta_roughness;
    app.surf_beretta.normal_texture = &app.tex_beretta_normal;
    app.surf_beretta.metallic_texture = NULL;
    app.surf_beretta.emissive_texture = NULL;

    app.surf_light.albedo = glms_vec3_zero();
    app.surf_light.emissive = glms_vec3_scale( (vec3s){1.0f, 0.10f, 0.00f}, 10000.0f);
    app.surf_light.ior = 0.0f;
    app.surf_light.metallic = 0.0f;
    app.surf_light.roughness= 0.0f;
    app.surf_light.material_func= material_light;
    app.surf_light.albedo_texture = NULL;
    app.surf_light.roughness_texture = NULL;
    app.surf_light.metallic_texture = NULL;
    app.surf_light.emissive_texture = NULL;

    float roughness_values[5] = {
        0.04f,
        0.2f,
        0.3f,
        0.4f,
        0.5f
    };
    float roughness_power = 1.0f;
    for(uint32_t i =0 ; i < 5; i++)
    {
        app.surf_gold[i].material_func = material_default;
        app.surf_gold[i].metallic = 1.0f;
        app.surf_gold[i].albedo = (vec3s){0.944f,0.776f,0.373f};
        app.surf_gold[i].emissive = glms_vec3_zero();
        app.surf_gold[i].ior = 0.47f;
        app.surf_gold[i].roughness = powf(roughness_values[i], roughness_power);
        app.surf_gold[i].albedo_texture = NULL;
        app.surf_gold[i].roughness_texture = NULL;
        app.surf_gold[i].metallic_texture = NULL;
        app.surf_gold[i].emissive_texture = NULL;

        app.surf_plastic[i].material_func = material_default;
        app.surf_plastic[i].metallic = 0.0f;
        app.surf_plastic[i].albedo = (vec3s){0.97f, 0.93f, 0.87f};
        app.surf_plastic[i].emissive = glms_vec3_zero();
        app.surf_plastic[i].ior = 1.45f;
        app.surf_plastic[i].roughness = powf(roughness_values[i], roughness_power);
        app.surf_plastic[i].albedo_texture = NULL;
        app.surf_plastic[i].roughness_texture = NULL;
        app.surf_plastic[i].metallic_texture = NULL;
        app.surf_plastic[i].emissive_texture = NULL;
    }
}

int32_t main(int32_t argc, char **argv)
{
    app.render_size = (ivec2s){1280, 720};
    app.fps_frames = 0;
    app.fps_seconds = 0.0;

    sw_initialize(app.render_size.x, app.render_size.y, 1.0f);

    setup_assets();



    // create camera
    camera_init(&app.camera, glm_rad(65.0f), (float)app.render_size.x / (float)app.render_size.y, 0.01f, 10.0f);
    app.pitch =  glm_rad(0.0f);
    versors pitch_rotation = glms_quatv(app.pitch, world_right);
    app.yaw = glm_rad(0.0f);
    versors yaw_rotation = glms_quatv(app.yaw, world_up);
    app.camera.rotation = glms_quat_mul(yaw_rotation, pitch_rotation);
    app.camera.position.z = 2.0f;
    app.camera.whitepoint = 20.0f;
    camera_apply(&app.camera);

    // create scene
    scene_create(&app.scene);
    app.scene.environment_hdr = &app.cube_env;
    app.scene.environment_intensity = 10.0f;

    entity e;

    e.mesh = &app.mesh_beretta;
    e.surface = &app.surf_beretta;
    entity_set_transform(&e, (vec3s){0.0f, 0.0f, 0.0f},  glms_quatv(glm_rad(90.0f), world_right));
    array_push(&app.scene.entities, &e);


    float step_size = glm_rad(360.0f / 5.0f);
    float radius = 1.0f;
    for(uint32_t i = 0; i < 5; i++)
    {
        float cx = cosf(i * step_size) * radius;
        float cy = sinf(i * step_size) * radius;

        cx = i * 1.5f;
        cy = 0.0f;

        e.mesh = &app.mesh_bust;
        e.surface = &app.surf_gold[i];
        entity_set_transform(&e, (vec3s){cx, 0.0f, -2.0f + cy},  glms_quatv(glm_rad(90.0f), world_right));
        array_push(&app.scene.entities, &e);

        e.mesh = &app.mesh_bust;
        e.surface = &app.surf_plastic[i];
        entity_set_transform(&e, (vec3s){cx, 0.0f, -4.0f + cy}, glms_quatv(glm_rad(90.0f), world_right));
        array_push(&app.scene.entities, &e);
    }

   
    app.tracer_cfg.camera = &app.camera; 
    app.tracer_cfg.scene = &app.scene;
    app.tracer_cfg.framebuffer_size = app.render_size;

    app.tracer_cfg.tile_size = (ivec2s){32, 32};
    app.tracer_cfg.num_threads = 15;
    //tracer_cfg_preset_high(&app.tracer_cfg);
    app.tracer_cfg.min_reflected_rays = 7;
    app.tracer_cfg.max_reflected_rays = 12;
    app.tracer_cfg.num_diffuse_rays = 16;
    app.tracer_cfg.max_bounces = 2;
    app.tracer_cfg.num_samples = 12;

    tracer_create(&app.tracer, &app.tracer_cfg);

    while(!sw_wants_quit())
    {
        sw_tick();
        double delta_time = sw_delta_time();

        app.seconds_since_preview += delta_time;

        bool is_rendering = app.tracer.rendering;
        bool was_rendering = is_rendering;
        if(is_rendering && tracer_render_is_done(&app.tracer))
        {
            is_rendering = false;
            tracer_render_end(&app.tracer);
        }

        if(!is_rendering)
        {
            float move_speed = 2.0f;
            float move_delta = delta_time * move_speed;

            float view_speed = 1.0f;
            float view_delta = delta_time * view_speed;

            if(sw_key_down(key_w))
            {
                app.camera.position = glms_vec3_add(app.camera.position, glms_vec3_scale(app.camera.derivatives.forward, move_delta));
            }

            if(sw_key_down(key_s))
            {
                app.camera.position = glms_vec3_add(app.camera.position, glms_vec3_scale(app.camera.derivatives.forward, -move_delta));
            }

            if(sw_key_down(key_d))
            {
                app.camera.position = glms_vec3_add(app.camera.position, glms_vec3_scale(app.camera.derivatives.right, move_delta));
            }

            if(sw_key_down(key_a))
            {
                app.camera.position = glms_vec3_add(app.camera.position, glms_vec3_scale(app.camera.derivatives.right, -move_delta));
            }

            if(sw_key_down(key_space))
            {
                app.camera.position = glms_vec3_add(app.camera.position, glms_vec3_scale(world_up, move_delta));
            }

            if(sw_key_down(key_lctrl))
            {
                app.camera.position = glms_vec3_add(app.camera.position, glms_vec3_scale(world_up, -move_delta));
            }

            if(sw_key_down(key_up))
            {
                app.pitch -= view_delta;
            }

            if(sw_key_down(key_down))
            {
                app.pitch += view_delta;
            }

            if(sw_key_down(key_right))
            {
                app.yaw += view_delta;
            }
            
            if(sw_key_down(key_left))
            {
                app.yaw -= view_delta;
            }

            versors pitch_rotation = glms_quatv(app.pitch, world_right);
            versors yaw_rotation = glms_quatv(app.yaw, world_up);
            app.camera.rotation = glms_quat_mul(yaw_rotation, pitch_rotation);

            camera_apply(&app.camera);
        }

        if(sw_key_pressed(key_r))
        {
            if(!is_rendering)
            {
                printf("begin rendering\n");
                tracer_render_begin(&app.tracer);
            }
        }

        if(sw_key_pressed(key_escape))
        {
            if(is_rendering)
            {
                printf("cancel rendering\n");
                tracer_render_cancel(&app.tracer);
            }
        }

        // calculate and print fps every second
        app.fps_seconds += delta_time;
        if(app.fps_seconds > 1.0)
        {
            printf("fps: %u\n", app.fps_frames);
            app.fps_frames = 0;
            app.fps_seconds = 0.0;
        }
        app.fps_frames++;



        if(is_rendering)
        {
            color *win_fb = sw_begin_render();
            vec3s *trace_fb = tracer_lock_framebuffer(&app.tracer);
            for(uint32_t y = 0; y < app.render_size.y; y++)
            {
                for(uint32_t x = 0; x < app.render_size.x; x++)
                {
                    uint32_t index = x + y * app.render_size.x;
                    win_fb[index] = vec3_to_color(trace_fb[index]);
                }
            }
            tracer_unlock_framebuffer(&app.tracer);
            sw_end_render();

            color indicator_color = (color){255, 255, 255, 220};
            color indicator_color2 = (color){0, 0, 0, 220};

            color indicator_color3 = (color){255, 255, 255, 80};
    
            for(int32_t thread_id = 0; thread_id < app.tracer.cfg->num_threads; thread_id++)
            {
                tracer_thread_indicator *indicator = &app.tracer.indicators[thread_id];
                if(!indicator->active)
                    continue;

                vec2s offset = {indicator->offset.x + 3.0f, indicator->offset.y + 3.0f};
                vec2s size = {indicator->size.x - 6.0f, indicator->size.y- 6.0f};
            
                sw_render_quad(offset, size,indicator_color3);



                vec2s tl = offset;
                vec2s tr = glms_vec2_add(offset, (vec2s){size.x, 0.0f});
                vec2s br = glms_vec2_add(offset, (vec2s){size.x, size.y});
                vec2s bl = glms_vec2_add(offset, (vec2s){0.0f, size.y});

                sw_render_line(tl, tr, indicator_color2);
                sw_render_line(tl, bl, indicator_color2);
                sw_render_line(tr, br, indicator_color2);
                sw_render_line(bl, br, indicator_color2);


                offset = (vec2s){indicator->offset.x + 4.0f, indicator->offset.y + 4.0f};
                size = (vec2s){indicator->size.x - 8.0f, indicator->size.y- 8.0f};

                tl = offset;
                tr = glms_vec2_add(offset, (vec2s){size.x, 0.0f});
                br = glms_vec2_add(offset, (vec2s){size.x, size.y});
                bl = glms_vec2_add(offset, (vec2s){0.0f, size.y});

                sw_render_line(tl, tr, indicator_color);
                sw_render_line(tl, bl, indicator_color);
                sw_render_line(tr, br, indicator_color);
                sw_render_line(bl, br, indicator_color);

            }

            struct timespec time;
            time.tv_sec = 0;
            time.tv_nsec = 100000000;


            thrd_sleep(&time, NULL);
        }
        else 
        {
            color *win_fb = sw_begin_render();

            if(was_rendering)
            {
                vec3s *trace_fb = tracer_lock_framebuffer(&app.tracer);
                for(uint32_t y = 0; y < app.render_size.y; y++)
                {
                    for(uint32_t x = 0; x < app.render_size.x; x++)
                    {
                        uint32_t index = x + y * app.render_size.x;
                        win_fb[index] = vec3_to_color(trace_fb[index]);
                    }
                }
                tracer_unlock_framebuffer(&app.tracer);
            }

            render_viewport(win_fb);
            sw_end_render();
        }

        sw_present();
    }


    if(app.tracer.rendering)
        tracer_render_cancel(&app.tracer);

    texture_destroy(&app.cube_env);
    mesh_destroy(&app.mesh_bust);
    mesh_destroy(&app.mesh_cube);
    mesh_destroy(&app.mesh_beretta);
    texture_destroy(&app.tex_beretta_albedo);
    texture_destroy(&app.tex_beretta_normal);
    texture_destroy(&app.tex_beretta_roughness);
    scene_destroy(&app.scene);
    tracer_destroy(&app.tracer);
    sw_shutdown();


}