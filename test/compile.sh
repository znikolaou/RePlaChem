#!/bin/bash

SRC_UTIL=$WRK_DIR/char_util.f90 
SRC_CHEM=$WRK_DIR/parse_chem.f90
SRC_GLBL=$WRK_DIR/Global.f90 $WRK_DIR/prec.f90
echo $SRC_CHEM

# $FC $FOPT $SRC_GLBL testIsBolsigReaction.f90 $SRC_CHEM -o goTestIsBolsig
# $FC $FOPT $SRC_GLBL testGetReactionSpecies.f90 $SRC_CHEM -o goTestGetReactionSpecies
# $FC $FOPT $SRC_GLBL testFormatReactions.f90 $SRC_CHEM -o goTestFormatReactions
# $FC $FOPT $SRC_GLBL testSeparateReaction.f90 $SRC_CHEM -o goTestSeparate
# $FC $FOPT $SRC_GLBL testCountSpeciesReactantsProducts.f90 $SRC_CHEM-o goTestCount

