#!/bin/bash
if [ -z "$1" ]
  then
    echo "Give name of folder to be searhced for images to format"
fi

#find $1 -type f -exec echo "put {}" \;
files=($(find $1 -type f))
for file in "${files[@]}"
do
    #echo "$file"
    if [ "${file: -4}" == ".jpg" ] || [ "${file: -4}" == ".png" ]
    then
        echo "$file : Converting image to ppm"
        mogrify -format ppm -quality 100 $file #2> /dev/null
        rm "$file"
    elif [ "${file: -4}" == ".mtl" ]
    then
        if grep -q -e jpg -e png "$file"
        then
            echo "$file : Changing image references to ppm"
            sed -i -e 's/jpg/ppm/g' -e 's/png/ppm/g' "$file"
        fi
    fi
done

echo "All files in $1 formatted correctly"