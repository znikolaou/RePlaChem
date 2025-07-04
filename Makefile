#
# AUTHOR: Z. NIKOLAOU
#

SRC_IN=global.f90 chem_parse.f90 io.f90\
       graph_search.f90 main.f90  \
       drg.f90 char_util.f90 timing.f90 \
       sort.f90 set.f90 stats.f90

SRC=$(addprefix $(REDCHEM_SRC),$(notdir $(SRC_IN)))

OBJS=$(patsubst $(REDCHEM_SRC)%.o,$(REDCHEM_BUILD)%.o,$(SRC:%.f90=%.o)) 

rePlaChem:$(OBJS) 
	$(FC) $(FOPT) $(OBJS) -o $(REDCHEM_BIN)/$@

$(REDCHEM_BUILD)%.o: $(REDCHEM_SRC)%.f90
	$(FC) $(FOPT) -c $< -o $@

clean:
	@ rm -f $(REDCHEM_BIN)/* *.o *.mod $(REDCHEM_BUILD)/*.o 

printObjs:
	@ echo $(OBJS)

printSrc:
	@ echo $(SRC)


