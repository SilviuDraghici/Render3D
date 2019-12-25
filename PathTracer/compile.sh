#!/bin/sh

# Use this command to compile single-threaded for developing and debugging your code
#g++ -O3 -g svdDynamic.c PathTracer.c utils_path.c -lm -o PathTracer

# Use this command to compile multi-threaded so you can quickly render complex images
# once your code is working
g++ -O3 -g -fopenmp svdDynamic.c PathTracer.c utils_path.c -lm -o PathTracer
