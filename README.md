# swgpu - A software path tracer written in C11 

* Cook Torrance Specular with GGX microfacet distribution
* Custom Diffuse distribution using the Golden Spiral, inspired by how sunflowers distribute their seeds evenly in the flower (converges much quicker than random hemisphere noise)
* Educational project, not for production use
* Only implements BRDF (specular, diffuse), not BSDF (transmittance, refraction etc)
* Supports custom materials, includes GoldenSpiral/CookTorrance BRDF and emissive materials
* Supports FBX meshes
* Supports any image formats that stb_image support, and supports bilinear texture sampling (no mipmaps yet)
