
cmake_minimum_required(VERSION 3.24)

project(swgpu)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${CMAKE_SOURCE_DIR}/bin)

set(SDL_SHARED OFF CACHE BOOL "Build Shared SDL3" FORCE)
set(SDL_STATIC ON CACHE BOOL "Build Static SDL3" FORCE)
add_subdirectory("ext/SDL" EXCLUDE_FROM_ALL)

set(CGLM_SHARED OFF CACHE BOOL "Build Shared CGLM" FORCE)
set(CGLM_STATIC ON CACHE BOOL "Build Static CGLM" FORCE)
add_subdirectory("ext/cglm" EXCLUDE_FROM_ALL)

add_library(swgpu.core STATIC "core/core.h" "core/core.c" "core/utils.h" "core/utils.c")
set_property(TARGET swgpu.core PROPERTY C_STANDARD 11)

target_include_directories(swgpu.core PUBLIC "core/")
if(SDL_SHARED)
    list(APPEND SWGPU_PRIVATE_DEPS SDL3::SDL3-shared)
elseif(SDL_STATIC)
    list(APPEND SWGPU_PRIVATE_DEPS SDL3::SDL3-static)
endif() 

list(APPEND SWGPU_PUBLIC_DEPS cglm)

target_link_libraries(swgpu.core PUBLIC ${SWGPU_PUBLIC_DEPS} PRIVATE ${SWGPU_PRIVATE_DEPS})


add_executable(swgpu.raytracer "raytracer/main.c"
    "raytracer/bxdf.c" "raytracer/bxdf.h"
    "raytracer/camera.c" "raytracer/camera.h"
    "raytracer/texture.c" "raytracer/texture.h"
    "raytracer/mesh.c" "raytracer/mesh.h"
    "raytracer/scene.c" "raytracer/scene.h"
    "raytracer/common.c" "raytracer/common.h"
    "raytracer/pathtracer.c" "raytracer/pathtracer.h"
    "raytracer/ufbx.h" "raytracer/ufbx.c"
    "raytracer/stb_image.h" "raytracer/stb_image.c"
    "raytracer/mtwister.h" "raytracer/mtwister.c")
set_property(TARGET swgpu.raytracer PROPERTY C_STANDARD 11)
target_link_libraries(swgpu.raytracer swgpu.core)