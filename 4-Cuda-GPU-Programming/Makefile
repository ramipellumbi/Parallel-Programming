# This Makefile assumes the following module files are loaded:
#
# CUDA
#
# This Makefile will only work if executed on a GPU node.
#

NVCC = nvcc
NVCCFLAGS = -O3 -g -I include
LFLAGS = -lm -g

# Compiler-specific flags (by default, we always use sm_75)
GENCODE_SM75 = -gencode=arch=compute_75,code=sm_75
GENCODE = $(GENCODE_SM75)

SRC_DIR 	= src
OBJ_DIR 	= obj
BIN_DIR 	= bin

.SUFFIXES : .cu .o

CU_FILES = $(wildcard $(SRC_DIR)/*.cu)
C_FILES = $(wildcard $(SRC_DIR)/*.c)
OBJ_FILES_CU = $(patsubst $(SRC_DIR)/%.cu, $(OBJ_DIR)/%.o, $(CU_FILES))
OBJ_FILES_C = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(C_FILES))
OBJ_FILES = $(OBJ_FILES_CU) $(OBJ_FILES_C)
BINARIES_CU = $(patsubst $(SRC_DIR)/%.cu, $(BIN_DIR)/%, $(CU_FILES))
BINARIES_C = $(patsubst $(SRC_DIR)/%.c, $(BIN_DIR)/%, $(C_FILES))
BINARIES = $(BINARIES_CU) $(BINARIES_C)


all: $(BINARIES)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	@mkdir -p $(@D)
	$(NVCC) $(GENCODE) $(NVCCFLAGS) -o $@ -c $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

$(BIN_DIR)/%: $(OBJ_DIR)/%.o
	@echo "Building binary $@ from object $^"
	@mkdir -p $(BIN_DIR)
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $^

%: $(BIN_DIR)/% 
	@:

# Clean rule
clean:	
	rm -f $(OBJ_DIR)/*.o $(BINARIES)
