

FC = ifort
AR = xiar
FTARGET = fgensln

##All taken from Intel link line advisor 
FLININCLUDEFLAGS =  -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
LINMKLFLAGS =   ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
MACMKLFLAGS= ${MKLROOT}/lib/libmkl_lapack95_lp64.a -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
FMACINCLUDEFLAGS= -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
UNAME_S:= $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	FINCLUDEFLAGS = $(FLININCLUDEFLAGS)
	MKLFLAGS=$(LINMKLFLAGS)
endif
ifeq ($(UNAME_S),Darwin)
	FINCLUDEFLAGS=$(FMACINCLUDEFLAGS)
	MKLFLAGS=$(MACMKLFLAGS)
endif


SRCS_ORD = base_parameters combinatorialsubs modelparms blas_interface   zetas operatorsolutions sys_solve_mod c_interface

SRC_DIR = .
BUILD_DIR = fbuild
MOD_DIR = fmod
SRCS_NAMES = $(patsubst %,$(SRC_DIR)/%.f90,$(SRCS_ORD))
MODS_NAMES = $(patsubst %,$(MOD_DIR)/%.mod,$(SRCS_ORD))
OBJS_NAMES = $(patsubst %,$(BUILD_DIR)/%.o,$(SRCS_ORD))


FFLAGS := -fPIC -nofor-main -module $(MOD_DIR) $(FINCLUDEFLAGS)
#-check all -traceback -profile-loops=all 
all: FFLAGS := -warn all  -fPIC -xHost -ipo -O3 -fltconsistency  $(FINCLUDEFLAGS) -module $(MOD_DIR)
all: $(OBJS_NAMES)
	$(FC) $(FFLAGS) sys_solve_prog.f90 -o $(FTARGET) $(OBJS_NAMES) $(MKLFLAGS)
	xiar r fortlibs.a $(OBJS_NAMES)

debug: FFLAGS := -g -O0 $(FFLAGS)
debug: $(OBJS_NAMES)
	$(FC) $(FFLAGS) sys_solve_prog.f90 -o $(FTARGET) $(OBJS_NAMES) $(MKLFLAGS)
	xiar r fortlibs.a $(OBJS_NAMES)


opt: FFLAGS := -O3 -ipo $(FFLAGS)
opt: all

libsonly: $(OBJS_NAMES)
	@echo $(SRCS_NAMES)
	@echo "Building libs..."
	xiar -r fortlibs.a $(OBJS_NAMES)

# $(BUILD_DIR)/blas_interface.o: $(SRC_DIR)/blas_interface.f90
# 	$(FC) $(FFLAGS) $(FLAPACK95_INCLUDE) -c $^ -o $@

$(BUILD_DIR)/%.o $(MOD_DIR)/%.mod: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $^ -o $@

.PHONY: clean
clean:
	rm fmod/* fbuild/* $(FTARGET) fortlibs.a

