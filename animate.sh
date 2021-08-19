#!/bin/bash
NAME=${1:-out}
RES=256

echo_err() { printf "%s\n" "$*" >&2; }

echo_err Resolution: ${RES}x${RES}
echo_err

mkdir -p frames
rm -rf frames/*
#Recycle Bin
for i in `seq -w 1 239`
do
    echo_err "Frame $i"
    ./Render3D 0 $RES 100 zz.png $i
    mv zz.png frames/$i.png
    echo_err
done

#the stderr is piped to null so it doesn't get recorded when logging render information
mencoder "mf://frames/*.png" -mf type=png:fps=24 -ovc x264 -x264encopts subq=6:partitions=all:8x8dct:me=umh:frameref=5:bframes=3:b_pyramid=normal:weight_b -o $NAME.avi 2> /dev/null
ffmpeg -y -i $NAME.avi $NAME.gif 2> /dev/null