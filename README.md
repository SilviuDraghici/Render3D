# Render3D
![FC]

3D RayTracing Engine - A Ray Tacer and Path Tracer all-in-one
===
This program can render scenes using either Ray tracing and the Phong illumination model or Path tracing.

## Implemented

Geometry:
- Arbitrary object transformations
- Meshes from Wavefront OBJ files (.obj, .mtl)
- Basic L-Systems
- Primitive shapes: Spheres, Cubes, Planes, Triangles...

Shading:
- Ambient Occlusion
- Path tracing with next event estimation
- Direct Lighting

Materials (Path tracing):
- Specular reflection
- Specular transmission
- Diffuse

Acceleration:
- Bounding Volume Hierarchy (BVH)
- Explicit light sampling
- Cosine Weighted Sampling
- Russian roulette

Misc:
- Texture mapping
- Normal mapping
- Alpha mapping
- Color transform: Linear to sRGB

---
Compiling
---
Use the MakeFile to compile.
Just run `make`

Running
---
 Run `Render3D` to get instructions on what parameters need to be passed in.

TODO
---
- Use Cuda to render on the gpu
- Refactor code to look cleaner and use c++ features more consistently.
- Implement or use a math library with SSE instructions

Sample Renders
---

1024x1024 Path traced Cornell Box. 50,000 samples.
![PT] 

1024x1024 Ray traced Cornell Box.
![RT]

4096x4096 Path Traced Wineglass Mesh
![WG]

1024x1024 Video of Path Traced Wineglass Mesh
![WGG]

4096x4096 Pathtraced Scene of 1M Triangle Mesh \"Lucy\"
Meant to make use of BVH's
![LT]

---

Some meshes from the following resource have been tested successfully: (McGuire Computer Graphics Archive)[https://casual-effects.com/data/].


[FC]: https://github.com/SilviuDraghici/Render3D/raw/master/output/fancy_apartment/fancy_kitchen.png "This is my favorite render so far!"

[PT]: https://github.com/SilviuDraghici/Render3D/raw/master/output/Cornell_Box_50k.png "This looks no better than at 20,000 samples. ¯\\_( ツ )_/¯ "

[RT]: https://github.com/SilviuDraghici/Render3D/raw/master/output/Cornell_Box_rt.png "This is a little less feature-full"

[WG]: https://github.com/SilviuDraghici/Render3D/raw/master/output/wineglass.png

[WGG]: https://github.com/SilviuDraghici/Render3D/raw/master/output/wineglass.gif

[LT]: https://github.com/SilviuDraghici/Render3D/raw/master/output/Lucy_brings_flowers.png "She is holding a simple L-system generated flower bouquet"