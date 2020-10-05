3D Render - A Ray Tacer and Path Tracer all-in-one
===
This program can render scenes using either Ray tracing and the Phong illumination model or Path tracing.

The idea is that, since ray tracing is much faster to render, scenes or features can be tested with ray tracing and final renders can then be done with path tracing.

Compiling
---
Use the MakeFile to compile.
Just run `make`

Running
---
 Run `Render3D` to get instructions on what parameters need to be passed in.

TODO
---
- implement volumes
- Use Cuda to render on the gpu
- use Lua Scripting to handle buildscenes

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

4096x4096 Pathtraced Sceene of 1M triangle Mesh \"Lucy\"
Meant to make use of BVH's
![LT]

<!--- Images ---> 

[PT]: https://github.com/SilviuDraghici/Render3D/raw/master/output/Cornell_Box_50k.png "This looks no better than at 20,000 samples. ¯\\_( ツ )_/¯ "

[RT]: https://github.com/SilviuDraghici/Render3D/raw/master/output/Cornell_Box_rt.png "This is a little less feature-full"

[WG]: https://github.com/SilviuDraghici/Render3D/raw/master/output/wineglass.png

[WGG]: https://github.com/SilviuDraghici/Render3D/raw/master/output/wineglass.gif

[LT]: https://github.com/SilviuDraghici/Render3D/raw/master/output/Lucy_brings_flowers.png "She is holding a simple L-system generated flower bouquet"