#!/bin/sh

####################
# Program to generate test images to test our segmentation
# script, and then run the segmentation process on those test
# images
####################

#Constants
EXPECTED_ARGS=4
E_BADARGS=65
DIRECTORY=test_segmentation

#check command line args
if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: $0 {im_dim} {shape_dim} {offset} {grid_size}"
    exit $E_BADARGS
fi

#Print command line args
echo "Test image dimensions: $1"
echo "Test shape dimensions: $2"
echo "Test image offset: $3"
echo "Test image grid size: $4"

#Make test image directory if doesnt exist
if [ ! -d "$DIRECTORY" ]; then
    mkdir $DIRECTORY
fi


echo "Generating test images"

#Generate basic ellipsoid image
vgenim x=$1 y=$1 z=$1 e=$2,$2,$2 | vdim -c of=$DIRECTORY/in_circle_normal_nogrid.vx

#Offset basic ellipsoid
vgenim x=$1 y=$1 z=$1 e=$2,$2,$2 l=$3,$3,$3 | vdim -c of=$DIRECTORY/in_circle_offset_nogrid.vx

#Basic ellipsoid with grid
vgenim x=$1 y=$1 z=$1 e=$2,$2,$2 xg=$4 yg=$4 zg=$4 | vdim -c of=$DIRECTORY/in_circle_normal_grid.vx

#Offset basic ellipsoid with grid
vgenim x=$1 y=$1 z=$1 e=$2,$2,$2 l=$3,$3,$3 xg=$4 yg=$4 zg=$4 | vdim -c of=$DIRECTORY/in_circle_offset_grid.vx

#Basic rectangle image
vgenim x=$1 y=$1 z=$1 r=$2,$2,$2 | vdim -c of=$DIRECTORY/in_rect_normal_nogrid.vx

#Offset basic rectangle
vgenim x=$1 y=$1 z=$1 r=$2,$2,$2 l=$3,$3,$3 | vdim -c of=$DIRECTORY/in_rect_offset_nogrid.vx

#Basic rectangle image with grid
vgenim x=$1 y=$1 z=$1 r=$2,$2,$2 xg=$4 yg=$4 zg=$4 | vdim -c of=$DIRECTORY/in_rect_normal_grid.vx

#Offset basic rectangle with grid
vgenim x=$1 y=$1 z=$1 r=$2,$2,$2 l=$3,$3,$3 xg=$4 yg=$4 zg=$4 | vdim -c of=$DIRECTORY/in_rect_offset_grid.vx
echo

echo "Running segmentation script on tests"

./segment.sh $DIRECTORY/in_circle_normal_nogrid.vx 2 1 $DIRECTORY/out_circle_normal_nogrid.vx

./segment.sh $DIRECTORY/in_circle_offset_nogrid.vx 2 1 $DIRECTORY/out_circle_offset_nogrid.vx

./segment.sh $DIRECTORY/in_circle_normal_grid.vx 2 1 $DIRECTORY/out_circle_normal_grid.vx

./segment.sh $DIRECTORY/in_circle_offset_grid.vx 2 1 $DIRECTORY/out_circle_offset_grid.vx

./segment.sh $DIRECTORY/in_rect_normal_nogrid.vx 2 1 $DIRECTORY/out_rect_normal_nogrid.vx

./segment.sh $DIRECTORY/in_rect_offset_nogrid.vx 2 1 $DIRECTORY/out_rect_offset_nogrid.vx

./segment.sh $DIRECTORY/in_rect_normal_grid.vx 2 1 $DIRECTORY/out_rect_normal_grid.vx

./segment.sh $DIRECTORY/in_rect_offset_grid.vx 2 1 $DIRECTORY/out_rect_offset_grid.vx
echo

echo "Calculating Volumes"

echo "Sphere Volumes"
./v3dvol if=$DIRECTORY/out_circle_normal_nogrid.vx
./v3dvol if=$DIRECTORY/out_circle_offset_nogrid.vx
./v3dvol if=$DIRECTORY/out_circle_normal_grid.vx
./v3dvol if=$DIRECTORY/out_circle_offset_grid.vx
echo "Cube Volumes"
./v3dvol if=$DIRECTORY/out_rect_normal_nogrid.vx
./v3dvol if=$DIRECTORY/out_rect_offset_nogrid.vx
./v3dvol if=$DIRECTORY/out_rect_normal_grid.vx
./v3dvol if=$DIRECTORY/out_rect_offset_grid.vx

