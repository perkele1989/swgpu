
#include "core.h"

#define SDL_MAIN_HANDLED
#include <SDL3/SDL.h>

typedef struct sw_ctx
{
    bool initialized;
    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_Texture* framebuffer;
    uint64_t ns;
    uint64_t delta_ns;
    double delta_time;
    bool wants_quit;

    uint64_t frame_index;

    uint16_t width;
    uint16_t height;
    double render_scaling;
    uint16_t render_width;
    uint16_t render_height;

    bool keys_state[key_count];
    bool keys_pressed[key_count];
    bool keys_released[key_count];
} sw_ctx;

static sw_ctx g_ctx = { 0 };


vec4s color_to_vec4(color in)
{
    vec4s c = {in.r, in.g, in.b, in.a};
    return glms_vec4_scale(c, 1.0f / 255.0f);   
}

color vec4_to_color(vec4s in)
{
    vec4s v = glms_vec4_scale(glms_vec4_clamp(in, 0.0f, 1.0f), 255.0f);
    return (color){(uint8_t)v.r, (uint8_t)v.g, (uint8_t)v.b, (uint8_t)v.a};
}

color vec3_to_color(vec3s in)
{
    vec3s v = glms_vec3_scale(glms_vec3_clamp(in, 0.0f, 1.0f), 255.0f);
    return (color){(uint8_t)v.r, (uint8_t)v.g, (uint8_t)v.b, (uint8_t)255};
}


sw_status sw_initialize(uint16_t w, uint16_t h, double render_scaling)
{
    if(g_ctx.initialized)
    {
        return sw_already_initialized;
    }

    g_ctx.frame_index = 0;
    g_ctx.width = w;
    g_ctx.height = h;
    g_ctx.render_scaling = render_scaling;
    g_ctx.render_width = w * render_scaling;
    g_ctx.render_height = h * render_scaling;

    if(!SDL_InitSubSystem(SDL_INIT_VIDEO|SDL_INIT_EVENTS))
    {
        printf("SDL_InitSubSystem failed: %s\n", SDL_GetError());
        return sw_fail;
    }

    if(!SDL_CreateWindowAndRenderer("swgpu", g_ctx.width, g_ctx.height, 0, &g_ctx.window, &g_ctx.renderer))
    {
        printf("SDL_CreateWindowAndRenderer failed: %s\n", SDL_GetError());
        SDL_QuitSubSystem(SDL_INIT_VIDEO|SDL_INIT_EVENTS);
        return sw_fail;
    }

    g_ctx.framebuffer = SDL_CreateTexture(g_ctx.renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STREAMING, g_ctx.render_width, g_ctx.render_height);
    if(!g_ctx.framebuffer)
    {
        printf("SDL_CreateTexture for framebuffer failed: %s", SDL_GetError());
        SDL_DestroyRenderer(g_ctx.renderer);
        SDL_DestroyWindow(g_ctx.window);
        SDL_QuitSubSystem(SDL_INIT_VIDEO|SDL_INIT_EVENTS);
        return sw_fail;
    }

    SDL_SetRenderDrawBlendMode(g_ctx.renderer, SDL_BLENDMODE_BLEND);
    SDL_SetTextureScaleMode(g_ctx.framebuffer, SDL_SCALEMODE_LINEAR);

    g_ctx.initialized = true;

    return sw_ok;
}

void sw_tick()
{
    bool want_quit = false;
    SDL_Event event;
    while (SDL_PollEvent(&event))
    {
        if (event.type == SDL_EVENT_QUIT)
        {
            want_quit = true;
        }
    }

    if(want_quit)
    {
        g_ctx.wants_quit = want_quit;
    }

    uint64_t ns = SDL_GetTicksNS();
    g_ctx.delta_ns = ns - g_ctx.ns;
    g_ctx.ns = ns;
    g_ctx.delta_time = ((double)g_ctx.delta_ns) / 1000000000.0;
    g_ctx.frame_index++;
    bool *sdl_key_state = SDL_GetKeyboardState(NULL);

    for(uint32_t i = 0; i < key_count; i++)
    {
        bool old_state = g_ctx.keys_state[i];
        bool new_state = sdl_key_state[i];
        g_ctx.keys_state[i] = new_state;
        g_ctx.keys_pressed[i] = new_state && !old_state;
        g_ctx.keys_released[i] = !new_state && old_state;
    }

}
ivec2s sw_render_size()
{
    return (ivec2s){ .x = g_ctx.render_width, .y = g_ctx.render_height};
}



color *sw_begin_render()
{
    void* pixels = NULL;
    int32_t pitch = 0;
    if(!SDL_LockTexture(g_ctx.framebuffer, NULL, &pixels, &pitch ))
    {
        printf("SDL_LockTexture failed: %s\n", SDL_GetError());
        return NULL;
    }

    if(pitch != g_ctx.render_width*4)
    {
        printf("Pitch mismatch\n");
    }

    return pixels;
}

void sw_end_render()
{
    SDL_UnlockTexture(g_ctx.framebuffer);

    SDL_SetRenderDrawColor(g_ctx.renderer, 255, 255, 255, 255);
    SDL_RenderTexture(g_ctx.renderer, g_ctx.framebuffer, NULL, NULL);
}

void sw_present()
{
    SDL_RenderPresent(g_ctx.renderer);
}

void sw_skip_render()
{
    SDL_SetRenderDrawColor(g_ctx.renderer, 255, 255, 255, 255);
    SDL_RenderTexture(g_ctx.renderer, g_ctx.framebuffer, NULL, NULL);
}


void sw_render_quad(vec2s offset, vec2s size, color col)
{
    SDL_SetRenderDrawColor(g_ctx.renderer, col.r, col.g, col.b, col.a);
    SDL_FRect rect;
    rect.w = size.x;
    rect.h = size.y;
    rect.x = offset.x;
    rect.y = offset.y;
    SDL_RenderFillRect(g_ctx.renderer, &rect);
}

void sw_render_line(vec2s a, vec2s b, color col)
{
    SDL_SetRenderDrawColor(g_ctx.renderer, col.r, col.g, col.b, col.a);
    SDL_RenderLine(g_ctx.renderer, a.x, a.y, b.x, b.y);
}

void sw_shutdown()
{
    if(!g_ctx.initialized)
        return;
    SDL_Quit();
    g_ctx.initialized = true;
}


bool sw_wants_quit()
{
    bool returner = g_ctx.wants_quit;
    g_ctx.wants_quit = false;
    return returner;
}

double sw_delta_time()
{
    return g_ctx.delta_time;
}


double sw_seconds()
{
    return ((double)g_ctx.ns) / 1000000000.0;
}

bool sw_key_down(key k)
{
    return g_ctx.keys_state[k];
}

bool sw_key_pressed(key k)
{
    return g_ctx.keys_pressed[k];
}

bool sw_key_released(key k)
{
    return g_ctx.keys_released[k];
}

uint64_t sw_frame_index()
{
    return g_ctx.frame_index;
}