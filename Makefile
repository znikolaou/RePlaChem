
SRC=$(SRC_DIR)/global.f90 $(SRC_DIR)/io.f90 $(SRC_DIR)/chem_parse.f90 \
    $(SRC_DIR)/graph_search.f90 $(SRC_DIR)/main.f90  \
    $(SRC_DIR)/drg.f90 $(SRC_DIR)/output.f90 $(SRC_DIR)/charUtil.f90 $(SRC_DIR)/timing.f90 \
    $(SRC_DIR)/sort.f90

OBJS=$(patsubst $(SRC_DIR)/%.o,$(BUILD_DIR)/%.o,$(SRC:%.f90=%.o)) 

redChem:$(OBJS) 
	$(FC) $(FOPT) $(OBJS) -o $(BIN_DIR)/$@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FOPT) -c $< -o $@

clean:
	@ rm -f $(BIN_DIR)/* *.o *.mod $(BUILD_DIR)/*.o 

print:
	@ echo $(OBJS)


