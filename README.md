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
- Re-implement meshes
- refactor with more c++ features
- Use Vulkan to render on the gpu

Sample Renders
---

1024x1024 Path traced Cornell Box. 50,000 samples.
[![PT](https://github.com/SilviuDraghici/Render3D/tree/master/output/Cornell_Box_50k.png)](This looks no better than at 20,000 samples. ¯\\_( ツ )_/¯ )

1024x1024 Ray traced Cornell Box. 
[![PT](https://github.com/SilviuDraghici/Render3D/tree/master/output/Cornell_Box_rt.png)](This looks no better than at 20,000 samples. ¯\\_( ツ )_/¯ )


