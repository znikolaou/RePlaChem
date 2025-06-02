
SRC_CHEM=global.f90 chem_parse.f90 char_util.f90 \
	 io.f90 timing.f90 sort.f90 
SRC=$(addprefix $(REDCHEM_SRC),$(notdir $(SRC_CHEM)))
OBJS_CHEM=$(patsubst $(REDCHEM_SRC)%.f90,$(REDCHEM_BUILD)%.o,$(SRC))

SRC_T1=testChemParse.f90

OBJS_T1=$(OBJS_CHEM) $(SRC_T1:.f90=.o)

printSourceFiles: 	
	@echo $(SRC_T1)

printObjectFiles: 	
	@echo $(OBJS_T1)

testChemParse: $(OBJS_T1)
	$(FC) $(FOPT) $(OBJS_T1) -o $@


%.o: %.f90
	$(FC) $(FOPT) -c -o $@ $< 

$(REDCHEM_BUILD)%.o: $(REDCHEM_SRC)%.f90
	$(FC) $(FOPT) -c -o $@ $< 

clean:
	@ rm -f testChemParse $(REDCHEM_BUILD)*.o  *.o log* go*  *.mod

