#!/bin/bash
# Author: Z. Nikolaou (2025)
# Shell script for setting env. vars.
#

# Set below Fortran compiler and compiler options:
FC='gfortran'
FOPT='-O3 -mcmodel=large -fbounds-check'

#----------------------------------------------------------------------#
#--------------------No editing below this line -----------------------#

WRK_DIR=$(pwd)
SRC_DIR=$WRK_DIR/src/
BUILD_DIR=$WRK_DIR/build/
BIN_DIR=$WRK_DIR/bin/

echo " "
echo " **********plasRedChem_v0.1**********"
echo " Plasma-chemistry reduction library. "
echo " Author: Z. Nikolaou (2025)"
echo " ************************************"
echo " "
echo " Setting environment vars ..."
echo " Work dir:"
echo " "$WRK_DIR
echo " Build dir:"
echo " "$BUILD_DIR
echo " Bin dir:"
echo " "$BIN_DIR
echo " Fortran compiler:"
echo " "$FC
echo " Compiler options:"
echo " "$FOPT
echo ""
echo " Environment variables set!"
echo " "
echo " Please make clean and make redChem to compile."
echo " "

export PATH=$PATH:$BIN_DIR
export WRK_DIR
export SRC_DIR
export FC
export FOPT
export BUILD_DIR
export BIN_DIR
