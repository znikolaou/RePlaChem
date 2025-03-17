#!/bin/bash

WRK_DIR=$(pwd)
SRC_DIR=$WRK_DIR/src/
BUILD_DIR=$WRK_DIR/build/

FC='gfortran'
FOPT='-O3 -mcmodel=large -fbounds-check'

echo " "
echo "Setting environment vars"
echo " "
echo "Work dir:"
echo $WRK_DIR
echo "Fortran compiler:"
echo $FC
echo "Compiler options:"
echo $FOPT
echo "Build dir:"
echo $BUILD_DIR

export WRK_DIR
export SRC_DIR
export FC
export FOPT
export BUILD_DIR
