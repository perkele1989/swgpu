
#include "pathtracer.h"
#include "texture.h"
#include <string.h>

void tracer_cfg_preset_production(tracer_cfg *out_cfg)
{
    out_cfg->max_reflected_rays = 32;
    out_cfg->min_reflected_rays = 128;
    out_cfg->num_diffuse_rays = 256;
    out_cfg->max_bounces = 8;
    out_cfg->num_samples = 128;
}

void tracer_cfg_preset_high(tracer_cfg *out_cfg)
{
    out_cfg->max_reflected_rays = 16;
    out_cfg->min_reflected_rays = 4;
    out_cfg->num_diffuse_rays = 64;
    out_cfg->max_bounces = 4;
    out_cfg->num_samples = 32;
}

void tracer_cfg_preset_preview(tracer_cfg *out_cfg)
{
    out_cfg->max_reflected_rays = 7;
    out_cfg->min_reflected_rays = 3;
    out_cfg->num_diffuse_rays = 16;
    out_cfg->max_bounces = 2;
    out_cfg->num_samples = 8;
}

#define max_tile_size 256

vec3s tonemap_reinhard2(vec3s x, float w)
{
    const float ls = w*w;
    vec3s xoverls = glms_vec3_adds(glms_vec3_divs(x, ls), 1.0f);
  
    return glms_vec3_div(glms_vec3_mul(x, xoverls), glms_vec3_adds(x, 1.0f));
}


int32_t _tracer_thread_func(void* ctx)
{
    static thread_local vec3s tile_buffer[max_tile_size*max_tile_size];

    tracer_thread *thread = ctx;
    int32_t w = thread->tracer->cfg->framebuffer_size.x;
    int32_t h = thread->tracer->cfg->framebuffer_size.y;
    camera *cam = thread->tracer->cfg->camera;

    tracer_job job;
    while(tracer_fetch_job(thread->tracer, thread, &job))
    {
        for(int32_t sample_index = 0; sample_index < thread->tracer->cfg->num_samples; sample_index++)
        {
            thread->sample_index = sample_index;
            for(int32_t y = job.offset.y; y  < job.offset.y + job.size.y; y++)
            {
                for(int32_t x = job.offset.x; x < job.offset.x + job.size.x; x++)
                {
                    if(thread->kill)
                    {
                        goto on_kill;
                    }
                    int32_t tx = x - job.offset.x;
                    int32_t ty = y - job.offset.y;
                    int32_t tile_index = tx + ty * job.size.x;

                    float u = (float)x / (float)w;
                    float v = (float)y / (float)h;
                
                    float ndc_x = u * 2.0f - 1.0f;
                    float ndc_y = v * 2.0f - 1.0f;
                
                    vec3s view_vector = camera_view_vector(cam, ndc_x, ndc_y);
                
                    ray initial_ray;
                    initial_ray.direction = view_vector;
                    initial_ray.origin = cam->position;
                    
                    vec3s output = tracer_thread_trace(thread, initial_ray, thread->tracer->cfg->max_bounces + 1);
                
                    output = tonemap_reinhard2(output, cam->whitepoint);

                    if(sample_index == 0)
                    {
                        // first sample sets initial average
                        tile_buffer[tile_index] = output;
                    }
                    else 
                    {
                        // all other samples average to buffer
                        tile_buffer[tile_index] = glms_vec3_divs( glms_vec3_add(glms_vec3_scale(tile_buffer[tile_index], sample_index), output) , sample_index + 1.0f);
                    }
                }
            }

            // output to main framebuffer 
            mtx_lock(&thread->tracer->framebuffer_lock);
            for(int32_t y = job.offset.y; y  < job.offset.y + job.size.y; y++)
            {
                for(int32_t x = job.offset.x; x < job.offset.x + job.size.x; x++)
                {
                    int32_t tx = x - job.offset.x;
                    int32_t ty = y - job.offset.y;
                    int32_t tile_index = tx + ty * job.size.x;
                    int32_t framebuffer_index = x + y * w;
                    thread->tracer->framebuffer[framebuffer_index] = tile_buffer[tile_index];
                }
            }
            mtx_unlock(&thread->tracer->framebuffer_lock);
        }
    }

on_kill:
    thread->done = true;
    
    printf("killing thread %u\n", thread->index);

    return 0;
}

int32_t _cmp_tile_offset(void* ctx, void const* lhs, void const* rhs)
{
    ivec2s *center = ctx;
    ivec2s *left = lhs;
    ivec2s *right = rhs;

    vec2s lf = {left->x, left->y};
    vec2s rf = {right->x, right->y};
    vec2s centerf = {center->x, center->y};

    float len_l = glms_vec2_norm(glms_vec2_sub(centerf, lf));
    float len_r =  glms_vec2_norm(glms_vec2_sub(centerf, rf));
    
    if (len_l < len_r)
        return -1;
    else if (len_l > len_r)
        return 1;
    
    return 0;
}

void interleave_in_place(ivec2s* arr, int32_t num_jobs, int32_t num_threads)
{

    ivec2s* temp = malloc(num_jobs * sizeof(ivec2s));

    int32_t index = 0;
    for (int32_t t = 0; t < num_threads; t++) {
        for (int32_t i = t; i < num_jobs; i += num_threads) {
            temp[index++] = arr[i];
        }
    }

    memcpy(arr, temp, num_jobs * sizeof(ivec2s));
    free(temp);
}

void tracer_create(tracer *tracer, tracer_cfg *cfg)
{
    tracer->cfg = cfg;
    tracer->rendering = false;
    // framebuffer
    mtx_init(&tracer->framebuffer_lock, mtx_plain);
    tracer->framebuffer = malloc(sizeof(vec3s) * tracer->cfg->framebuffer_size.x * tracer->cfg->framebuffer_size.y);

    int32_t num_tiles_x = ceilf((float)tracer->cfg->framebuffer_size.x / (float)tracer->cfg->tile_size.x);
    int32_t num_tiles_y = ceilf((float)tracer->cfg->framebuffer_size.y / (float)tracer->cfg->tile_size.y);

    int32_t last_tile_width = tracer->cfg->framebuffer_size.x % tracer->cfg->tile_size.x;
    if (last_tile_width == 0)
        last_tile_width = tracer->cfg->tile_size.x;
    
        int32_t last_tile_height = tracer->cfg->framebuffer_size.y % tracer->cfg->tile_size.y;
    if (last_tile_height == 0)
        last_tile_height = tracer->cfg->tile_size.y;
    
    ivec2 last_tile_size = { last_tile_width, last_tile_height };

    int32_t num_tiles = num_tiles_x * num_tiles_y;

    ivec2s center_tile = { num_tiles_x/2, num_tiles_y/2 };

    int32_t tiles_per_thread = num_tiles / tracer->cfg->num_threads;
    int32_t threads_with_extra = num_tiles % tracer->cfg->num_threads;

    ivec2s *tile_indices = malloc(num_tiles * sizeof(ivec2s));
    for(int32_t y = 0; y < num_tiles_y; y++)
    {
        for(int32_t x = 0; x < num_tiles_x; x++)
        {
            int32_t array_index = y * num_tiles_x + x;
            tile_indices[array_index] = (ivec2s){x, y};
        }
    }

    qsort_s(tile_indices, num_tiles, sizeof(ivec2s), _cmp_tile_offset, &center_tile);

    // interleave it for better locality
    // interleave_in_place(tile_indices, num_tiles, cfg->num_threads);

    tracer->num_jobs = num_tiles;
    tracer->jobs = malloc(sizeof(tracer_job) * num_tiles);

    for(int32_t tile_index = 0; tile_index < num_tiles; tile_index++)
    {
        ivec2s tile = tile_indices[tile_index];

        tracer_job new_job;
        new_job.offset = glms_ivec2_mul(tile, tracer->cfg->tile_size);
        new_job.size = tracer->cfg->tile_size;

        if(tile.x == num_tiles_x - 1)
        {
            new_job.size.x = last_tile_width;
        }

        if(tile.y == num_tiles_y - 1)
        {
            new_job.size.y = last_tile_height;
        }

        tracer->jobs[tile_index] = new_job;
    }

    // then when distributing: tiles for that thread is =  tiles_per_thread + (thread_index < threads_with_extra ? 1 : 0)

    // create thread states and distribute jobs (we won't create the actual threads, thats done in tracer_render_begin)
    tracer->threads = malloc(sizeof(tracer_thread) * cfg->num_threads);
    tracer->indicators = malloc(sizeof(tracer_thread_indicator) * cfg->num_threads);
    for(uint32_t thread_index = 0; thread_index < cfg->num_threads; thread_index++)
    {
        tracer_thread *thread = &tracer->threads[thread_index];
        tracer_thread_indicator * indicator = &tracer->indicators[thread_index];
        indicator->active = false;

        thread->tracer = tracer;
        thread->index = thread_index;
        thread->rand = seedRand(13.042 + 0.252 * 323.3);
        array_create(&thread->path_stack, sizeof(path_state), 32);
        array_reserve(&thread->path_stack, tracer->cfg->max_bounces * (tracer->cfg->max_reflected_rays + tracer->cfg->num_diffuse_rays));

        thread->sample_index = 0;
        thread->done = false;
        thread->kill = false;
    }

    free(tile_indices);
}


bool tracer_fetch_job(tracer* tracer, tracer_thread *thread, tracer_job *out_job)
{
    mtx_lock(&tracer->jobs_lock);

    if(tracer->next_free_job >= tracer->num_jobs)
    {
        tracer->indicators[thread->index].active = false;
        mtx_unlock(&tracer->jobs_lock);
        return false;
    }

    memcpy(out_job, tracer->jobs + tracer->next_free_job, sizeof(tracer_job));
    tracer->next_free_job++;

    tracer->indicators[thread->index].active = true;
    tracer->indicators[thread->index].offset = out_job->offset;
    tracer->indicators[thread->index].size = out_job->size;

    mtx_unlock(&tracer->jobs_lock);

    return true;
}


void tracer_render_begin(tracer *tracer)
{
    if(tracer->rendering)
    {
        printf("already rendering\n");
        return;
    }
    tracer->rendering = true;

    // clear framebuffer 
    for(uint32_t y = 0; y < tracer->cfg->framebuffer_size.y; y++)
    {
        for(uint32_t x = 0; x < tracer->cfg->framebuffer_size.x; x++)
        {
            uint32_t index = x + y * tracer->cfg->framebuffer_size.x;
            tracer->framebuffer[index] = glms_vec3_zero();
        }
    }

    tracer->next_free_job = 0;

    // setup all threads 
    for(uint32_t thread_index = 0; thread_index < tracer->cfg->num_threads; thread_index++)
    {
        tracer_thread *thread = &tracer->threads[thread_index];
        tracer_thread_indicator *indicator = &tracer->indicators[thread_index];
        indicator->active = false;
        thread->done = false;
        thread->kill = false;
        thrd_create(&thread->handle, _tracer_thread_func, thread);
    }
}

bool tracer_render_is_done(tracer *tracer)
{
    if(!tracer->rendering)
    {
        return true;
    }

    for(uint32_t thread_index = 0; thread_index < tracer->cfg->num_threads; thread_index++)
    {
        tracer_thread *thread = &tracer->threads[thread_index];
        if(!thread->done)
        {
            return false;
        }
    }

    return true;
}

vec3s *tracer_render_end(tracer* tracer)
{
    if(!tracer->rendering)
    {
        printf("not rendering, will return NULL\n");
        return NULL;
    }

    tracer->rendering = false;
    tracer->next_free_job = 0;
    for(uint32_t thread_index = 0; thread_index < tracer->cfg->num_threads; thread_index++)
    {
        tracer_thread *thread = &tracer->threads[thread_index];
        tracer_thread_indicator *indicator = &tracer->indicators[thread_index];
        indicator->active = false;
        thrd_join(thread->handle, NULL);
        thread->kill = false;
        thread->done = false;
    }

    return tracer->framebuffer;
}

void tracer_render_cancel(tracer* tracer)
{
    if(!tracer->rendering)
    {
        printf("not rendering\n");
        return;
    }

    tracer->rendering = false;
    tracer->next_free_job = 0;
    for(uint32_t thread_index = 0; thread_index < tracer->cfg->num_threads; thread_index++)
    {
        tracer_thread *thread = &tracer->threads[thread_index];
        thread->kill = true;
        thrd_join(thread->handle, NULL);
    }
}


vec3s *tracer_lock_framebuffer(tracer *tracer)
{
    mtx_lock(&tracer->framebuffer_lock);
    return tracer->framebuffer;
}

void tracer_unlock_framebuffer(tracer *tracer)
{
    mtx_unlock(&tracer->framebuffer_lock);
}

void tracer_destroy(tracer *tracer)
{
    for(uint32_t thread_index = 0; thread_index < tracer->cfg->num_threads; thread_index++)
    {
        tracer_thread *thread = &tracer->threads[thread_index];
        if(tracer->rendering)
        {
            thread->kill = true;
            thrd_join(thread->handle, NULL);
        }
        array_destroy(&thread->path_stack);
    }

    free(tracer->indicators);
    free(tracer->jobs);

    free(tracer->framebuffer);
    free(tracer->threads);
}


vec3s tracer_thread_trace(tracer_thread *thread, ray input_ray, int32_t max_depth)
{
    // clear the path stack without reallocating
    array_clear(&thread->path_stack);

    // push the input ray to the path stack with a weight of one
    array_push(&thread->path_stack, &(path_state){input_ray, {1.0f, 1.0f, 1.0f}, max_depth});

    // initialize accumulator to zero
    vec3s accumulator = {0.0f, 0.0f, 0.0f};

    tracer_cfg* cfg = thread->tracer->cfg;

    // while we have paths in the stack ..
    while(thread->path_stack.size > 0)
    {
        if(thread->kill)
            return glms_vec3_zero();

        // pop the path at the top of the stack
        path_state *state = array_pop(&thread->path_stack);

        // we've reached our depth limit, so skip this ray
        if(state->depth <= 0)
            continue;

        // early out if weight is below treshold, avoids doing bounces that does not contribute, speeding things up
        if(glms_vec3_norm(state->weight) < 0.00001f)
            continue;

        // query the scene with the current ray (i.e.: shoot the ray)
        query query = scene_query(cfg->scene, state->ray);

        // no hit: sample environment
        if(!query.hit)
        {
            // sample the environment cubemap using ray direction as uvw
            vec3s hdr_color = glms_vec4_copy3(texture_sample_uvw(cfg->scene->environment_hdr, state->ray.direction)); 
            // hdr_color = srgb_to_linear(hdr_color);
            // add environment color weighted by intensity and path weight to accumulator
            accumulator = glms_vec3_add(accumulator, glms_vec3_mul(glms_vec3_scale(hdr_color, cfg->scene->environment_intensity), state->weight));
            continue;
        }

        // hit: invoke surface material function
        if(query.hit_surface->material_func)
        {
            // evaluate surface material
            vec3s emissive = query.hit_surface->material_func(thread, &query, state->weight, state->depth);

            // calculate inverse square attenuation (for emissive)
            float attenuation = 1.0f / (query.hit_distance * query.hit_distance + 1e-6f);

            // add emissive weighted by attentuation and path weight to accumulator
            accumulator = glms_vec3_add(accumulator, glms_vec3_scale(glms_vec3_mul(emissive, state->weight), attenuation));
        }
    }

    // return the accumulator
    return accumulator;
}


