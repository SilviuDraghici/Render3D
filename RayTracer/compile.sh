#!/bin/sh
g++ -O4 -g svdDynamic.c RayTracer.c utils.c affineTransforms.c -lm -o RayTracer
