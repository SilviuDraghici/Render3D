
# Render3D
![FK]

3D Ray Tracing Engine - A Ray Tacer and Path Tracer all-in-one
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
The program was built and tested on linux but it should work on mac and windows with little to no modification. 

Running
---
 Run `Render3D` in a shell to get complete instructions on what parameters need to be passed in.
 
 The first parameter is the rendering mode.
 - 0 - Ray tracing
 - 1 - Path tracing
 - 2 - normals debug
 - 3 - colors debug
<!-- -->
Example:

    ./Render3D 1 3840x2160 1000 out.png
This will render the current scene in path tracing mode at 3840x2160 resolution for 1000 samples and save it as out.png

TODO
---
- Use Cuda to ray trace on the gpu
- Change buildscene system to not require recompile for different scenes.
- Separate object color into diffuse and specular components.  (also ambient for phong illum)
- Refactor code to look cleaner and use c++ features more consistently.
- Implement or use a math library with SSE instructions

Sample Renders
---
3840x2160 Path traced living room scene. 1000 samples.
![LR_PT]

1920x1080 debug view displaying the object colors with no shading of any kind.
![LR_C]

1920x1080 debug view displaying the object normals.
![LR_N]

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

<!-- Images -->

[FK]: https://github.com/SilviuDraghici/Render3D/raw/master/output/fancy_apartment/fancy_kitchen.png "This is my favorite render so far!"

[PT]: https://github.com/SilviuDraghici/Render3D/raw/master/output/Cornell_Box_50k.png "This looks no better than at 20,000 samples. ¯\\_( ツ )_/¯ "

[RT]: https://github.com/SilviuDraghici/Render3D/raw/master/output/Cornell_Box_rt.png "This is a little less feature-full"

[WG]: https://github.com/SilviuDraghici/Render3D/raw/master/output/wineglass.png

[WGG]: https://github.com/SilviuDraghici/Render3D/raw/master/output/wineglass.gif "I had to turn this into a gif to get it to show up"

[LT]: https://github.com/SilviuDraghici/Render3D/raw/master/output/Lucy_brings_flowers.png "She is holding a simple L-system generated flower bouquet"

[LR_PT]: https://github.com/SilviuDraghici/Render3D/raw/master/output/living_room/living_room.png "I replaced the pictures on the wall with my own"

[LR_C]: https://github.com/SilviuDraghici/Render3D/raw/master/output/living_room/colors.png 

[LR_N]: https://github.com/SilviuDraghici/Render3D/raw/master/output/living_room/normals.png 
