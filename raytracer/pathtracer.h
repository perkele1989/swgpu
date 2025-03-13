
#pragma once 

#include "common.h"
#include "camera.h"
#include "scene.h"

#include "mtwister.h"


#include <core.h>
#include <threads.h>

/** configuration for the tracer */
typedef struct tracer_cfg
{
    // cant change after tracer_create
    ivec2s tile_size;

    // cant change after tracer_create
    ivec2s framebuffer_size;

    // cant change after tracer_create
    uint32_t num_threads;

    // camera used for rendering
    camera *camera;

    // scene used for rendering
    scene *scene;

    /** how many times a ray can bounce after initial ejection */
    uint32_t max_bounces;

    /** minimum number of reflection rays (at roughness = 0) */
    uint32_t min_reflected_rays;

    /** maximum number of reflection rays (at roughness = 1) */
    uint32_t max_reflected_rays;

    /** number of diffuse rays */
    uint32_t num_diffuse_rays;

    /** number of samples included in the final average */
    uint32_t num_samples;
    
} tracer_cfg;

void tracer_cfg_preset_production(tracer_cfg *out_cfg);
void tracer_cfg_preset_high(tracer_cfg *out_cfg);
void tracer_cfg_preset_preview(tracer_cfg *out_cfg);

/** state for one bounce in a traced path */
typedef struct path_state
{
    ray ray;
    vec3s weight;

    int32_t depth;
} path_state;

typedef struct tracer_job 
{
    ivec2s offset;
    ivec2s size;
} tracer_job;

typedef struct tracer_thread_indicator
{
    ivec2s offset;
    ivec2s size;
    bool active;
} tracer_thread_indicator;
 
typedef struct tracer_thread
{
    struct tracer *tracer;

    int32_t index;

    thrd_t handle;
    MTRand rand; 

    int32_t sample_index;

    bool kill;
    bool done;

    /** stack of path states for current fragment in current job (all the bounces) */
    array path_stack; // type: path_state
} tracer_thread;


vec3s tracer_thread_trace(tracer_thread *thread, ray input_ray, int32_t max_depth);

typedef struct tracer 
{
    tracer_cfg *cfg;

    // framebuffer lock, worker threads write into it, but main thread can read from it
    mtx_t framebuffer_lock;
    vec3s *framebuffer;
     

    tracer_job *jobs; 
    uint32_t num_jobs;
    uint32_t next_free_job;
    mtx_t jobs_lock;


    bool rendering;

    tracer_thread *threads;
    tracer_thread_indicator *indicators;
} tracer;

/** creates new tracer */
void tracer_create(tracer *tracer, tracer_cfg *cfg);

bool tracer_fetch_job(tracer* tracer, tracer_thread *thread, tracer_job *out_job);

/** begins rendering this tracer */
void tracer_render_begin(tracer *tracer);

/** returns true if tracer is done rendering */
bool tracer_render_is_done(tracer *tracer);

/** blocks execution until rendering is done, returns framebuffer */
vec3s *tracer_render_end(tracer* tracer);

/** cancel the render */
void tracer_render_cancel(tracer* tracer);

/** returns a locked framebuffer pointer, rendering may pause as long as this is locked */
vec3s *tracer_lock_framebuffer(tracer *tracer);

/** unlock a previously locked framebuffer */
void tracer_unlock_framebuffer(tracer *tracer);

/** destroys tracer */
void tracer_destroy(tracer *tracer);

vec3s tonemap_reinhard2(vec3s x, float w);