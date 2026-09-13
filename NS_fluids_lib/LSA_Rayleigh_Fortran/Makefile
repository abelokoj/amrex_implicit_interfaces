# =====================================================================
# LSA_Rayleigh_Fortran -- NS-MFP linear stability analysis
# =====================================================================
#   make            build the tools (bin/make_inputs, bin/run_lsa)
#   make test       build and run the whole test suite
#   make clean      remove build products
#
# Requires: a Fortran 2008 compiler and LAPACK/BLAS.
#   Debian/Ubuntu:  sudo apt-get install gfortran liblapack-dev libblas-dev
#   macOS:          brew install gcc  (LAPACK comes from the Accelerate
#                   framework; see LAPACK_LIBS below)
# ---------------------------------------------------------------------

# NOTE: FC is assigned with '=' not '?=' on purpose.  GNU make has a built-in
# default of FC=f77, and a built-in default BEATS a '?=' assignment, so '?='
# would silently select f77.  On Ubuntu f77 is a symlink to gfortran and it
# happens to work; on a cluster with a real f77 it fails confusingly.
# Override explicitly on the command line if needed:  make FC=ifx
FC      = gfortran
FFLAGS  ?= -O2 -g -Wall -Wno-maybe-uninitialized -std=f2008 -fimplicit-none

# LAPACK/BLAS.  On macOS use:  make LAPACK_LIBS="-framework Accelerate"
LAPACK_LIBS ?= -llapack -lblas

SRC   := src
TST   := tests
BUILD := build
BIN   := bin
MODS  := $(BUILD)/mods

# Module compilation order matters in Fortran: a module must be compiled
# before anything that USEs it.
LIB_SRC := \
	$(SRC)/table_mod.f90 \
	$(SRC)/dispersion_mod.f90 \
	$(SRC)/plotfile_mod.f90 \
	$(SRC)/dmd_mod.f90 \
	$(SRC)/nsmfp_mod.f90 \
	$(SRC)/proportionality_mod.f90 \
	$(SRC)/interface_mode_mod.f90 \
	$(SRC)/bundle_mod.f90 \
	$(SRC)/make_inputs_mod.f90

LIB_OBJ := $(patsubst $(SRC)/%.f90,$(BUILD)/%.o,$(LIB_SRC))

TEST_LIB_SRC := $(TST)/testing_mod.f90 $(TST)/plotfile_writer_mod.f90
TEST_LIB_OBJ := $(patsubst $(TST)/%.f90,$(BUILD)/%.o,$(TEST_LIB_SRC))

PROGS := $(BIN)/make_inputs $(BIN)/run_lsa $(BIN)/dump_amplitude \
         $(BIN)/dump_profiles
TESTS := $(BIN)/test_dispersion $(BIN)/test_dmd \
         $(BIN)/test_proportionality $(BIN)/test_bundle \
         $(BIN)/test_end_to_end

.PHONY: all test clean dirs

all: dirs $(PROGS)

dirs:
	@mkdir -p $(BUILD) $(MODS) $(BIN)

# --- library objects --------------------------------------------------
$(BUILD)/%.o: $(SRC)/%.f90 | dirs
	$(FC) $(FFLAGS) -J$(MODS) -I$(MODS) -c $< -o $@

$(BUILD)/%.o: $(TST)/%.f90 | dirs
	$(FC) $(FFLAGS) -J$(MODS) -I$(MODS) -c $< -o $@

# Explicit ordering: gfortran needs each .mod present before its users.
$(BUILD)/plotfile_mod.o:      $(BUILD)/dispersion_mod.o
$(BUILD)/run_lsa.o:           $(BUILD)/table_mod.o
$(BUILD)/make_inputs.o:       $(BUILD)/table_mod.o
$(BUILD)/dump_amplitude.o:    $(BUILD)/table_mod.o
$(BUILD)/dmd_mod.o:           $(BUILD)/dispersion_mod.o $(BUILD)/table_mod.o
$(BUILD)/nsmfp_mod.o:         $(BUILD)/plotfile_mod.o
$(BUILD)/proportionality_mod.o: $(BUILD)/dispersion_mod.o $(BUILD)/table_mod.o
$(BUILD)/interface_mode_mod.o: $(BUILD)/plotfile_mod.o
$(BUILD)/bundle_mod.o:        $(BUILD)/plotfile_mod.o
$(BUILD)/dump_profiles.o:     $(BUILD)/bundle_mod.o $(BUILD)/interface_mode_mod.o
$(BUILD)/make_inputs_mod.o:   $(BUILD)/dispersion_mod.o
$(BUILD)/plotfile_writer_mod.o: $(BUILD)/dispersion_mod.o

# --- programs ---------------------------------------------------------
$(BIN)/make_inputs: $(LIB_OBJ) $(SRC)/make_inputs.f90 | dirs
	$(FC) $(FFLAGS) -J$(MODS) -I$(MODS) $^ -o $@ $(LAPACK_LIBS)

$(BIN)/run_lsa: $(LIB_OBJ) $(SRC)/run_lsa.f90 | dirs
	$(FC) $(FFLAGS) -J$(MODS) -I$(MODS) $^ -o $@ $(LAPACK_LIBS)

$(BIN)/dump_amplitude: $(LIB_OBJ) $(SRC)/dump_amplitude.f90 | dirs
	$(FC) $(FFLAGS) -J$(MODS) -I$(MODS) $^ -o $@ $(LAPACK_LIBS)

$(BIN)/dump_profiles: $(LIB_OBJ) $(SRC)/dump_profiles.f90 | dirs
	$(FC) $(FFLAGS) -J$(MODS) -I$(MODS) $^ -o $@ $(LAPACK_LIBS)

# --- tests ------------------------------------------------------------
$(BIN)/test_%: $(LIB_OBJ) $(TEST_LIB_OBJ) $(TST)/test_%.f90 | dirs
	$(FC) $(FFLAGS) -J$(MODS) -I$(MODS) $^ -o $@ $(LAPACK_LIBS)

# test_end_to_end and test_proportionality write files under /tmp and read
# inputs/inputs_rayleigh_template.txt, so they must run from the project root.
test: dirs $(TESTS)
	@echo "======================================================"
	@echo " running the test suite"
	@echo "======================================================"
	@fail=0; \
	for t in $(TESTS); do \
	  echo ""; echo ">>> $$t"; \
	  ./$$t || fail=1; \
	done; \
	echo ""; \
	if [ $$fail -eq 0 ]; then \
	  echo "ALL TEST PROGRAMS PASSED"; \
	else \
	  echo "SOME TESTS FAILED"; exit 1; \
	fi

clean:
	rm -rf $(BUILD) $(BIN)
