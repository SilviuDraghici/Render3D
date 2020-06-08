#!/bin/bash
NAME=${1:-out}

mkdir -p frames
rm -rf frames/*
#Recycle Bin
for i in `seq -w 1 240`
do
    echo "Frame $i"
    ./Render3D 0 1024 100 out.ppm $i
    mv out.ppm frames/$i.ppm
    echo
done

mogrify -format jpg -quality 100 frames/*.ppm
rm frames/*.ppm
mencoder "mf://frames/*.jpg" -mf type=jpeg:fps=24 -ovc x264 -x264encopts subq=6:partitions=all:8x8dct:me=umh:frameref=5:bframes=3:b_pyramid=normal:weight_b -o $NAME.avi
ffmpeg -y -i $NAME.avi $NAME.gif