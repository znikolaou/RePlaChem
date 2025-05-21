
SRC_IN=global.f90 io.f90 chem_parse.f90 \
       graph_search.f90 main.f90  \
       drg.f90 output.f90 charUtil.f90 timing.f90 \
       sort.f90 set.f90

SRC=$(addprefix $(SRC_DIR),$(notdir $(SRC_IN)))

OBJS=$(patsubst $(SRC_DIR)%.o,$(BUILD_DIR)%.o,$(SRC:%.f90=%.o)) 

redChem:$(OBJS) 
	$(FC) $(FOPT) $(OBJS) -o $(BIN_DIR)/$@

$(BUILD_DIR)%.o: $(SRC_DIR)%.f90
	$(FC) $(FOPT) -c $< -o $@

clean:
	@ rm -f $(BIN_DIR)/* *.o *.mod $(BUILD_DIR)/*.o 

printObjs:
	@ echo $(OBJS)

printSrc:
	@ echo $(SRC)


