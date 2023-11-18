CC_SERIAL = icc
CC_MPI = mpicc

CFLAGS     = -g -O3 -xHost -fno-alias -std=c99 -I include
CFLAGS_OMP = -g -O3 -xHost -fno-alias -std=c99 -qopenmp -I include  

SRC_DIR 	= src
OBJ_DIR 	= obj
BIN_DIR 	= bin
INCLUDE_DIR = include

MATMUL    = $(OBJ_DIR)/matmul.o
TIMING 	  = $(OBJ_DIR)/timing.o
UTILITIES = $(OBJ_DIR)/utilities.o
MODULES   = $(MATMUL) $(TIMING) $(UTILITIES)

EXECUTABLES = 1-serial

# For debugging
SRC_FILES = $(shell find $(SRC_DIR) -name '*.c')
OBJ_FILES = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILES))
$(info SRC_FILES = $(SRC_FILES))
$(info OBJ_FILES = $(OBJ_FILES))

all: $(addprefix $(BIN_DIR)/, $(EXECUTABLES))

# general file rule
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC_MPI) $(CFLAGS) -c $< -o $@

# Specific rule for serial
$(OBJ_DIR)/1-serial.o: $(SRC_DIR)/1-serial.c
	@mkdir -p $(@D)
	$(CC_SERIAL) $(CFLAGS) -c $< -o $@	

1-serial: $(BIN_DIR)/1-serial
$(BIN_DIR)/1-serial: $(OBJ_DIR)/1-serial.o $(MODULES)
	@mkdir -p $(BIN_DIR)
	$(CC_SERIAL) $(CFLAGS) $^ -o $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)