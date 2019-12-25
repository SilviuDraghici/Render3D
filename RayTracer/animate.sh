#!/bin/bash
mkdir -p anim
rm -rf anim/*
for i in `seq -w 1 120`
do
    ./RayTracer 512 2 0 out.ppm $i
    mv out.ppm anim/$i.ppm
done
mogrify -format jpg -quality 100 anim/*.ppm
rm anim/*.ppm
mencoder "mf://anim/*.jpg" -mf type=jpeg:fps=24 -ovc x264 -x264encopts subq=6:partitions=all:8x8dct:me=umh:frameref=5:bframes=3:b_pyramid=normal:weight_b -o ${1:output}.avi
ffmpeg -y -i ${1:output}.avi ${1:output}.gif