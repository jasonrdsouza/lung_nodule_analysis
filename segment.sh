#!/bin/sh

####################
#Program to perform segmentation and noise removal of
#3d images
####################

#Constants
EXPECTED_ARGS=4
E_BADARGS=65

#check command line args
if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: $0 {in} {k_size} {thresh} {out}"
    exit $E_BADARGS
fi

echo "Input file: $1"
echo "Kernel size: $2"
echo "Threshold value: $3"
echo "Output file: $4"

#convert from short to byte (vfix)
echo "converting to byte image"
vfix -byte if=$1 of=temp_byte

#convert to binary image
echo "converting to thresholded binary image"
vpix if=temp_byte th=$3 hi=255 lo=0 of=temp_binary

#morphological filtering (v3morph) "opening"
echo "performing morphological filtering"
v3morph if=temp_binary -ed t=s s=$2,$2,$2 of=temp_opened

#conversion to 3d format (vdim)
echo "converting to visionx 3d format"
vdim if=temp_opened -c of=$4

#remove temporary files
echo "removing temporary files"
#comment below line to examine the temporary files for debugging
rm temp_byte temp_binary temp_opened

echo
