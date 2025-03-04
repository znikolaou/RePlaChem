
#Source files
SRC=globalModule.f90 precModule.f90 graphSearchModule.f90 parsechem.f90 main.f90 readRates.f90 \
    drg.f90  output.f90 charUtil.f90 timing.f90 sort.f90

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

