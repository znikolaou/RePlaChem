#
# AUTHOR: Z. NIKOLAOU
#

SRC_CHEM=global.f90 chem_parse.f90 char_util.f90 \
	 io.f90 timing.f90 sort.f90 
SRC_UTIL=pcmo.f90

SRC=$(addprefix $(REDCHEM_SRC),$(notdir $(SRC_CHEM)))
OBJS_CHEM=$(patsubst $(REDCHEM_SRC)%.f90,$(REDCHEM_BUILD)%.o,$(SRC))

SRCU=$(addprefix $(REDCHEM_UTIL),$(notdir $(SRC_UTIL)))
OBJSU=$(patsubst $(REDCHEM_UTIL)%.f90,$(REDCHEM_BUILD)%.o,$(SRCU))

SRC_T1=testChemParse.f90
SRC_T2=testPCMO.f90

OBJS_T1=$(OBJS_CHEM) $(SRC_T1:.f90=.o)
OBJS_T2=$(OBJSU) $(SRC_T2:.f90=.o)

printSourceFiles: 	
	@echo $(SRC_T1)

printObjectFiles: 	
	@echo $(OBJS_T1)

testChemParse: $(OBJS_T1)
	$(FC) $(FOPT) $(OBJS_T1) -o $@

testPCMO: $(OBJS_T2)
	$(FC) $(FOPT) $(OBJS_T2) -o $@

%.o: %.f90
	$(FC) $(FOPT) -c -o $@ $< 

$(REDCHEM_BUILD)%.o: $(REDCHEM_SRC)%.f90 
	$(FC) $(FOPT) -c -o $@ $< 

$(REDCHEM_BUILD)%.o: $(REDCHEM_UTIL)%.f90 
	$(FC) $(FOPT) -c -o $@ $< 


clean:
	@ rm -f testChemParse testPCMO $(REDCHEM_BUILD)*.o  *.o log* go*  *.mod ./output/*

