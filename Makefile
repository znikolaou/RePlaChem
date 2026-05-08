#
# AUTHOR: Z. NIKOLAOU
#

SRC_IN=main.f90 io.f90 timing.f90 sort.f90 stats.f90

SRC=$(addprefix $(REDCHEM_SRC),$(notdir $(SRC_IN)))

OBJS=$(patsubst $(REDCHEM_SRC)%.o,$(REDCHEM_BUILD)%.o,$(SRC:%.f90=%.o)) 

rePlaChem:$(OBJS) 
	$(FC) $(FOPT) $(OBJS) $(REDCHEM_SRC)/char_util.a $(REDCHEM_SRC)/chmp_dr.a $(REDCHEM_SRC)/set.a -o $(REDCHEM_BIN)/$@


$(REDCHEM_BUILD)%.o: $(REDCHEM_SRC)%.f90
	$(FC) $(FOPT) -c $< -o $@

clean:
	@ rm -f $(REDCHEM_BIN)/* *.o *.mod $(REDCHEM_BUILD)/*.o 

printObjs:
	@ echo $(OBJS)

printSrc:
	@ echo $(SRC)


