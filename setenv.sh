#!/bin/bash
# Author: Z. Nikolaou (2025)
# Shell script for setting env. vars.
#

# Set below Fortran compiler and compiler options:
FC='gfortran'
FOPT='-O3 -mcmodel=large -fbounds-check'

#----------------------------------------------------------------------#
#--------------------No editing below this line -----------------------#

REDCHEM_HOME=$(pwd)
REDCHEM_SRC=$REDCHEM_HOME/src/
REDCHEM_BUILD=$REDCHEM_HOME/build/
REDCHEM_BIN=$REDCHEM_HOME/bin/

echo " "
echo " **********plasRedChem_v0.1**********"
echo " Plasma-chemistry reduction library. "
echo " Author: Z. Nikolaou (2025)"
echo " ************************************"
echo " "
echo " Setting environment vars ..."
echo " Home dir:"
echo " "$REDCHEM_HOME
echo " Source dir:"
echo " "$REDCHEM_HOME
echo " Build dir:"
echo " "$REDCHEM_BUILD
echo " Bin dir:"
echo " "$REDCHEM_BIN
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
export REDCHEM_HOME
export REDCHEM_SRC
export REDCHEM_BUILD
export REDCHEM_BIN
export FC
export FOPT
