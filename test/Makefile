
SRC_CHEM=$(REDCHEM_SRC)global.f90  \
	 $(REDCHEM_SRC)chem_parse.f90 \
	 $(REDCHEM_SRC)char_util.f90 \
	 $(REDCHEM_SRC)io.f90 $(REDCHEM_SRC)timing.f90 $(REDCHEM_SRC)sort.f90 

OBJS_CHEM=$(patsubst $(REDCHEM_SRC)%.f90,$(REDCHEM_BUILD)%.o,$(SRC_CHEM))

SRC_T1=testIsBolsigReaction.f90
SRC_T2=testFormatReactions.f90
SRC_T3=testSeparateReaction.f90
SRC_T4=testGetReactionSpeciesFromReaction.f90
SRC_T5=testCountSpeciesReactantsProducts.f90
SRC_T6=testReadChem.f90
SRC_T7=testIsChargedSpecies.f90
SRC_T8=testSplitString.f90
SRC_T9=testRemoveDuplicateReactions.f90
SRC_T10=testGetChargedSpecies.f90
SRC_T11=testGetSpeciesIndex.f90
SRC_T12=testChemParse.f90
SRC_T13=testReadLines.f90
SRC_T14=testSplitSpaceSeparated.f90
SRC_T15=testGetIndex.f90
SRC_T16=testExtractNu.f90

OBJS_T1=$(OBJS_CHEM) $(SRC_T1:.f90=.o)
OBJS_T2=$(OBJS_CHEM) $(SRC_T2:.f90=.o)
OBJS_T3=$(OBJS_CHEM) $(SRC_T3:.f90=.o)
OBJS_T4=$(OBJS_CHEM) $(SRC_T4:.f90=.o) 
OBJS_T5=$(OBJS_CHEM) $(SRC_T5:.f90=.o)
OBJS_T6=$(OBJS_CHEM) $(SRC_T6:.f90=.o)
OBJS_T7=$(OBJS_CHEM) $(SRC_T7:.f90=.o)
OBJS_T8=$(OBJS_CHEM) $(SRC_T8:.f90=.o) 
OBJS_T9=$(OBJS_CHEM) $(SRC_T9:.f90=.o)
OBJS_T10=$(OBJS_CHEM) $(SRC_T10:.f90=.o)
OBJS_T11=$(OBJS_CHEM) $(SRC_T11)
OBJS_T12=$(OBJS_CHEM) $(SRC_T12)
OBJS_T13=$(OBJS_CHEM) $(SRC_T13)
OBJS_T14=$(OBJS_CHEM) $(SRC_T14)
OBJS_T15=$(OBJS_CHEM) $(SRC_T15)
OBJS_T16=$(OBJS_CHEM) $(SRC_T16)

printSourceFiles: 	
	@echo $(SRC_CHEM)

printObjectFiles: 	
	@echo $(OBJS_CHEM)

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

goTestReadChem: $(OBJS_T6)
	$(FC) $(FOPT) $(OBJS_T6) -o $@

goTestIsChargedSpecies: $(OBJS_T7)
	$(FC) $(FOPT) $(OBJS_T7) -o $@

goTestSplitString: $(OBJS_T8)
	$(FC) $(FOPT) $(OBJS_T8) -o $@

goTestRemoveDuplicateReactions: $(OBJS_T9)
	$(FC) $(FOPT) $(OBJS_T9) -o $@

goTestGetChargedSpecies: $(OBJS_T10)
	$(FC) $(FOPT) $(OBJS_T10) -o $@

goTestGetSpecIndex: $(OBJS_T11)
	$(FC) $(FOPT) $(OBJS_T11) -o $@

goTestChemParse: $(OBJS_T12)
	$(FC) $(FOPT) $(OBJS_T12) -o $@

goTestReadLines: $(OBJS_T13)
	$(FC) $(FOPT) $(OBJS_T13) -o $@

goTestSplitSpaceSeparated: $(OBJS_T14)
	$(FC) $(FOPT) $(OBJS_T14) -o $@

goTestGetIndex: $(OBJS_T15)
	$(FC) $(FOPT) $(OBJS_T15) -o $@

goTestExtractNu: $(OBJS_T16)
	$(FC) $(FOPT) $(OBJS_T16) -o $@

print:
	@ echo $(OBJS_CHEM) 
%.o: %.f90
	$(FC) $(FOPT) -c -o $@ $< 

$(REDCHEM_BUILD)%.o: $(REDCHEM_SRC)%.f90
	$(FC) $(FOPT) -c -o $@ $< 

clean:
	@ rm -f $(REDCHEM_BUILD)*.o  *.o log* go*  *.mod

