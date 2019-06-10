NVCC=nvcc
TARGET=main
NVFLAGS=  -arch=sm_61 -x cu -I. -dc
NVDBGFLAGS= -Xcompiler -rdynamic -G -g
SRC_DIR=.
BUILD_DIR=build
CU_SRCS_NAMES= main.cu randomnumbers.cu util.cu parameters.c

OBJS:=$(patsubst %.cu,$(BUILD_DIR)/%.o,$(CU_SRCS_NAMES))
OBJS:=$(patsubst %.c,$(BUILD_DIR)/%.o,$(OBJS))


all: $(OBJS)
	$(NVCC) -arch=sm_61 $(OBJS) -o $(TARGET)
	


debug: NVFLAGS := -O0 $(NVDBGFLAGS) $(NVFLAGS)
debug: all	

info: NVFLAGS := --ptxas-options="--warn-on-spills --warn-on-local-memory-usage" $(NVFLAGS)
info: all

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVCC) $(NVFLAGS)  $^ -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(NVCC) $(NVFLAGS)  $^ -o $@


.PHONY: clean
clean:
	rm $(BUILD_DIR)/* $(TARGET)


