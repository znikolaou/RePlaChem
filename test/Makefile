
#Source files
SRC_CHEM=$(WRK_DIR)/Global.f90 $(WRK_DIR)/char_util.f90 $(WRK_DIR)/prec.f90 $(WRK_DIR)/parsechem.f90 
OBJS_CHEM=$(SRC_CHEM:%.f90=%.o)

SRC_T1=testIsBolsigReaction.f90
SRC_T2=testFormatReactions.f90
SRC_T3=testSeparateReaction.f90
SRC_T4=testGetReactionSpecies.f90
SRC_T5=testCountSpeciesReactantsProducts.f90
SRC_T6=testGetReactionSpeciesWithBolsig.f90
SRC_T7=testIsChargedSpecies.f90
SRC_T8=testSplitString.f90
SRC_T9=testRemoveDuplicateReactions.f90

OBJS_T1=$(SRC_T1:.f90=.o) $(OBJS_CHEM)
OBJS_T2=$(SRC_T2:.f90=.o) $(OBJS_CHEM)
OBJS_T3=$(SRC_T3:.f90=.o) $(OBJS_CHEM)
OBJS_T4=$(SRC_T4:.f90=.o) $(OBJS_CHEM)
OBJS_T5=$(SRC_T5:.f90=.o) $(OBJS_CHEM)
OBJS_T6=$(SRC_T6:.f90=.o) $(OBJS_CHEM)
OBJS_T7=$(SRC_T7:.f90=.o) $(OBJS_CHEM)
OBJS_T8=$(SRC_T8:.f90=.o) $(OBJS_CHEM)
OBJS_T9=$(SRC_T9:.f90=.o) $(OBJS_CHEM)

goTestIsBolsig: $(OBJS_T1)
	$(FC) $(FOPT) $(OBJS_T1) -o $@ 

goTestFormatReactions: $(OBJS_T2)
	$(FC) $(FOPT) $(OBJS_T2) -o $@

goTestSeparateReactions: $(OBJS_T3)
	$(FC) $(FOPT) $(OBJS_T3) -o $@

goTestGetReactionSpecies: $(OBJS_T4)
	$(FC) $(FOPT) $(OBJS_T4) -o $@

goTestCountSpeciesReactantsProducts: $(OBJS_T5)
	$(FC) $(FOPT) $(OBJS_T5) -o $@

goTestGetReactionSpeciesWithBolsig: $(OBJS_T6)
	$(FC) $(FOPT) $(OBJS_T6) -o $@

goTestIsChargedSpecies: $(OBJS_T7)
	$(FC) $(FOPT) $(OBJS_T7) -o $@

goTestSplitString: $(OBJS_T8)
	$(FC) $(FOPT) $(OBJS_T8) -o $@

goTestRemoveDuplicateReactions: $(OBJS_T9)
	$(FC) $(FOPT) $(OBJS_T9) -o $@

%.o: %.f90
	$(FC) $(FOPT) -c -o $@ $< 

$(WRK_DIR)/%.o: $(WRK_DIR)/%.f90
	$(FC) $(FOPT) -c -o $@ $< 

clean:
	rm $(WRK_DIR)/*.o  *.o log* go*

