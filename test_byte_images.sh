#!/bin/sh

####################
# Program to test our segmentation script on byte images
# as well as perform volume calculations
####################

#Constants
EXPECTED_ARGS=6
E_BADARGS=65
BIN_DIRECTORY=~/ece5470/project

#check command line args
if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: $0 {k_size} {bins} {img_dir} {xres} {yres} {zres}"
    exit $E_BADARGS
fi

DATA_DIRECTORY=$3
xres=$4
yres=$5
zres=$6

#go into test directory
cd $DATA_DIRECTORY
#clean out old results
rm *.out

echo "Running segmentation script on tests"
for file in `ls`; do
cd $BIN_DIRECTORY
./adaptive_segment_byte.sh $DATA_DIRECTORY/"$file" $1 $2 $DATA_DIRECTORY/"$file.out"
cd $DATA_DIRECTORY
done

echo "Begin Volume Calculation"
for file in `ls | grep .out`; do
echo "$file"
cd $BIN_DIRECTORY
./v3dvol_voxel if=$DATA_DIRECTORY/"$file" xres=$xres yres=$yres zres=$zres
cd $DATA_DIRECTORY
done

echo "End Volume Calculation"

