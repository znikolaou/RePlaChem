
#Source files
SRC=globalModule.f90 prec.f90 graph_search_mod.f90 parsechem.f90 main_redchem.f90 read_rates.f90 \
    drg.f90  output.f90 char_util.f90 timing.f90 sort.f90

OBJS=$(SRC:.f90=.o)

goRedChem: $(OBJS) 
	$(FC) $(FOPT) $(OBJS) -o $@

%.o: %.f90
	$(FC) $(FOPT) -c $< 


clean: clean_src clean_output

clean_src:
	rm exe* *.o *.mod log* goRedChem
clean_output:
	rm  output/*
	mkdir output 

