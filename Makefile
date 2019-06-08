CC=icc
TARGET=runsolver
BUILD_DIR=build
FORT_DIR=fort
#MKL FLAGS ARE FROM INTEL LINK LINE ADVISOR
UNAME_S:= $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	MKLFLAGS= -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm -ldl -lifcore ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a
endif
ifeq ($(UNAME_S),Darwin)
	MKLFLAGS= ${MKLROOT}/lib/libmkl_lapack95_lp64.a  -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/../compiler/lib/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lifcore
endif


SRCS_ORD = call_gensln generatevars f_interface
FORT_OBJS = $(wildcard fort/fbuild/*.o)

SRCS_NAMES = $(patsubst %,%.c,$(SRCS_ORD))
C_OBJS_NAMES = $(patsubst %,$(BUILD_DIR)/%.o,$(SRCS_ORD))

CFLAGS= -O3 -xHost -ipo

all: FORTLIBS
all: $(C_OBJS_NAMES)
	icc $(CFLAGS) -o $(TARGET) $(C_OBJS_NAMES)   $(FORT_OBJS)  $(MKLFLAGS)

debug: CFLAGS= -g -O0
debug: FORTLIBS_DEBUG
debug: $(C_OBJS_NAMES)
	icc $(CFLAGS) -o $(TARGET) $(C_OBJS_NAMES)  $(FORT_DIR)/fortlibs.a $(MKLFLAGS)

FORTLIBS:
	$(MAKE) -C $(FORT_DIR)

FORTLIBS_DEBUG:
	cd $(FORT_DIR) && $(MAKE) debug && cd ..


$(BUILD_DIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(MKLFLAGS)

.PHONY: clean
clean:
	rm build/* $(FORT_DIR)/fbuild/*.o $(FORT_DIR)/fmod/*.mod $(FORT_DIR)/fortlibs.a $(FORT_DIR)/fgensln $(TARGET)
