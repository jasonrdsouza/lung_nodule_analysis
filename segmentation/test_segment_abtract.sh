#!/bin/sh

####################
# Program to test our segmentation script on test images
# as well as perform volume calculations
####################

#Constants
EXPECTED_ARGS=2
E_BADARGS=65
DIRECTORY=abstract_images

#check command line args
if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: $0 {k_size} {thresh}"
    exit $E_BADARGS
fi


echo "Running segmentation script on tests"

./segment.sh $DIRECTORY/test1_cube.vx $1 $2 $DIRECTORY/out_test1_cube.vx

./segment.sh $DIRECTORY/test1_lin.vx $1 $2 $DIRECTORY/out_test1_lin.vx

./segment.sh $DIRECTORY/test1_nn.vx $1 $2 $DIRECTORY/out_test1_nn.vx

./segment.sh $DIRECTORY/test1_noisy.vx $1 $2 $DIRECTORY/out_test1_noisy.vx

./segment.sh $DIRECTORY/test1_noise_lin.vx $1 $2 $DIRECTORY/out_test1_noise_lin.vx

./segment.sh $DIRECTORY/test1_noise_nn.vx $1 $2 $DIRECTORY/out_test1_noise_nn.vx

./segment.sh $DIRECTORY/test2_sphere.vx $1 $2 $DIRECTORY/out_test2_sphere.vx

./segment.sh $DIRECTORY/test2_lin.vx $1 $2 $DIRECTORY/out_test2_lin.vx

./segment.sh $DIRECTORY/test2_nn.vx $1 $2 $DIRECTORY/out_test2_nn.vx

./segment.sh $DIRECTORY/test2_noisy.vx $1 $2 $DIRECTORY/out_test2_noisy.vx

./segment.sh $DIRECTORY/test2_noise_lin.vx $1 $2 $DIRECTORY/out_test2_noise_lin.vx

./segment.sh $DIRECTORY/test2_noise_nn.vx $1 $2 $DIRECTORY/out_test2_noise_nn.vx
echo


echo "Begin Volume Calculation"

echo "out_test1_cube.vx"
./v3dvol if=$DIRECTORY/out_test1_cube.vx
echo "out_test1_lin.vx"
./v3dvol if=$DIRECTORY/out_test1_lin.vx
echo "out_test1_nn.vx"
./v3dvol if=$DIRECTORY/out_test1_nn.vx
echo "out_test1_noisy.vx"
./v3dvol if=$DIRECTORY/out_test1_noisy.vx
echo "out_test1_noise_lin.vx"
./v3dvol if=$DIRECTORY/out_test1_noise_lin.vx
echo "out_test1_noise_nn.vx"
./v3dvol if=$DIRECTORY/out_test1_noise_nn.vx

echo "out_test2_sphere.vx"
./v3dvol if=$DIRECTORY/out_test2_sphere.vx
echo "out_test2_lin.vx"
./v3dvol if=$DIRECTORY/out_test2_lin.vx
echo "out_test2_nn.vx"
./v3dvol if=$DIRECTORY/out_test2_nn.vx
echo "out_test2_noisy.vx"
./v3dvol if=$DIRECTORY/out_test2_noisy.vx
echo "out_test2_noise_lin.vx"
./v3dvol if=$DIRECTORY/out_test2_noise_lin.vx
echo "out_test2_noise_nn.vx"
./v3dvol if=$DIRECTORY/out_test2_noise_nn.vx

echo "End Volume Calculation"

