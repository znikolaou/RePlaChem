
#Source files
SRC_CHEM=$(SRC_DIR)/globalModule.f90 $(SRC_DIR)/charUtil.f90 $(SRC_DIR)/precModule.f90 $(SRC_DIR)/parsechem.f90 
OBJS_CHEM=$(SRC_CHEM:%.f90=%.o)

SRC_T1=testIsBolsigReaction.f90
SRC_T2=testFormatReactions.f90
SRC_T3=testSeparateReaction.f90
SRC_T4=testGetReactionSpeciesFromReaction.f90
SRC_T5=testCountSpeciesReactantsProducts.f90
SRC_T6=testGetReactionSpeciesFromSpeciesListOnly.f90
SRC_T7=testIsChargedSpecies.f90
SRC_T8=testSplitString.f90
SRC_T9=testRemoveDuplicateReactions.f90
SRC_T10=testGetChargedSpecies.f90

OBJS_T1=$(SRC_T1:.f90=.o) $(OBJS_CHEM)
OBJS_T2=$(SRC_T2:.f90=.o) $(OBJS_CHEM)
OBJS_T3=$(SRC_T3:.f90=.o) $(OBJS_CHEM)
OBJS_T4=$(SRC_T4:.f90=.o) $(OBJS_CHEM)
OBJS_T5=$(SRC_T5:.f90=.o) $(OBJS_CHEM)
OBJS_T6=$(SRC_T6:.f90=.o) $(OBJS_CHEM)
OBJS_T7=$(SRC_T7:.f90=.o) $(OBJS_CHEM)
OBJS_T8=$(SRC_T8:.f90=.o) $(OBJS_CHEM)
OBJS_T9=$(SRC_T9:.f90=.o) $(OBJS_CHEM)
OBJS_T10=$(SRC_T10:.f90=.o) $(OBJS_CHEM)

goTestIsBolsig: $(OBJS_T1)
	$(FC) $(FOPT) $(OBJS_T1) -o $@ 

goTestFormatReactions: $(OBJS_T2)
	$(FC) $(FOPT) $(OBJS_T2) -o $@

goTestSeparateReactions: $(OBJS_T3)
	$(FC) $(FOPT) $(OBJS_T3) -o $@

goTestGetReactionSpeciesFromReaction: $(OBJS_T4)
	$(FC) $(FOPT) $(OBJS_T4) -o $@

goTestCountSpeciesReactantsProducts: $(OBJS_T5)
	$(FC) $(FOPT) $(OBJS_T5) -o $@

goTestGetReactionSpeciesFromSpeciesListOnly: $(OBJS_T6)
	$(FC) $(FOPT) $(OBJS_T6) -o $@

goTestIsChargedSpecies: $(OBJS_T7)
	$(FC) $(FOPT) $(OBJS_T7) -o $@

goTestSplitString: $(OBJS_T8)
	$(FC) $(FOPT) $(OBJS_T8) -o $@

goTestRemoveDuplicateReactions: $(OBJS_T9)
	$(FC) $(FOPT) $(OBJS_T9) -o $@

goTestGetChargedSpecies: $(OBJS_T10)
	$(FC) $(FOPT) $(OBJS_T10) -o $@

%.o: %.f90
	$(FC) $(FOPT) -c -o $@ $< 

$(SRC_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FOPT) -c -o $@ $< 

clean:
	rm $(SRC_DIR)/*.o  *.o log* go*  *.mod

