#!/bin/bash

$FC $FOPT $WRK_DIR/Global.f90 $WRK_DIR/prec.f90 testSeparateReaction.f90 $WRK_DIR/char_util.f90 $WRK_DIR/parsechem.f90 -o goTestSeparate
$FC $FOPT testCountSpeciesReactantsProducts.f90 $WRK_DIR/char_util.f90 -o goTestCount

