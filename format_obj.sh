#!/bin/bash
if [ -z "$1" ]
  then
    echo "Give name of .obj to be formatted"
fi

grep '# Vertices:' $1 > v2_$1
grep '# Faces:' $1 >> v2_$1
grep 'v ' $1 >> v2_$1
grep 'vn' $1 >> v2_$1
grep 'f ' $1 >> v2_$1
mv v2_$1 scenes/$1
rm $1