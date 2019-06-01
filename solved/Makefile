#A non-portable makefile for this project
#CC=icc
FC = ifort
# CTARGET = cgensln.so
FTARGET = fgensln
# CALLGEN = callgen
# CFLAGS = -I${MKLROOT}/include -fPIC
FINCLUDEFLAGS =  -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
MKLFLAGS =   ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

#to compile c from fort src
##something like icc -I/usr/include/libxml2 -o call_gensln.o call_gensln.c parameters.c my_lib.a ~/intelb/lib/intel64/libmkl_blas95_lp64.a ~/intelb/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core /opt/intel/compilers_and_libraries_2019.0.117/linux/compiler/lib/intel64_lin/libifcore.a -lxml2



SRCS_ORD =  modelparms blas_interface combinatorialsubs zetas operatorsolutions sys_solve_mod

SRC_DIR = .
BUILD_DIR = fbuild
MOD_DIR = fmod
SRCS_NAMES = $(patsubst %,$(SRC_DIR)/%.f90,$(SRCS_ORD))
MODS_NAMES = $(patsubst %,$(MOD_DIR)/%.mod,$(SRCS_ORD))
OBJS_NAMES = $(patsubst %,$(BUILD_DIR)/%.o,$(SRCS_ORD))


FFLAGS := -fPIC -nofor-main -module $(MOD_DIR) $(FINCLUDEFLAGS)

prog: FFLAGS := -fPIC -g -ipo -O2 -fltconsistency -profile-loops=all $(FINCLUDEFLAGS) -module $(MOD_DIR)
prog: $(OBJS_NAMES)
	$(FC) $(FFLAGS) sys_solve_prog.f90 -o $(FTARGET) $(OBJS_NAMES) $(MKLFLAGS)

debug: FFLAGS := -g -O0 $(FFLAGS)
debug: all

opt: FFLAGS := -O3 -ipo $(FFLAGS)
opt: all

all: $(OBJS_NAMES)
	@echo $(SRCS_NAMES)

$(BUILD_DIR)/blas_interface.o: $(SRC_DIR)/blas_interface.f90
	$(FC) $(FFLAGS) $(FLAPACK95_INCLUDE) $(FMKL_INCLUDE_NAMES) -c $^ -o $@

$(BUILD_DIR)/%.o $(MOD_DIR)/%.mod: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $^ -o $@

.PHONY: clean
clean:
	rm fmod/* fbuild/*

