CC_SERIAL = icc
CC_MPI = mpicc

CFLAGS = -g -O3 -xHost -fno-alias -std=c99 -I include

SRC_DIR 	= src
OBJ_DIR 	= obj
BIN_DIR 	= bin
INCLUDE_DIR = include

MATMUL    = $(OBJ_DIR)/matmul.o
TIMING 	  = $(OBJ_DIR)/timing.o
UTILITIES = $(OBJ_DIR)/utilities.o
MODULES   = $(MATMUL) $(TIMING) $(UTILITIES)

EXECUTABLES = serial task5 task6 task7 task8

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
$(OBJ_DIR)/serial.o: $(SRC_DIR)/serial.c
	@mkdir -p $(@D)
	$(CC_SERIAL) $(CFLAGS) -c $< -o $@	

serial: $(BIN_DIR)/serial
$(BIN_DIR)/serial: $(OBJ_DIR)/serial.o $(MODULES)
	@mkdir -p $(BIN_DIR)
	$(CC_SERIAL) $(CFLAGS) $^ -o $@

task5: $(BIN_DIR)/task5
$(BIN_DIR)/task5: $(OBJ_DIR)/task5.o $(MODULES)
	@mkdir -p $(BIN_DIR)
	$(CC_MPI) $(CFLAGS) $^ -o $@

task6: $(BIN_DIR)/task6
$(BIN_DIR)/task6: $(OBJ_DIR)/task6.o $(MODULES)
	@mkdir -p $(BIN_DIR)
	$(CC_MPI) $(CFLAGS) $^ -o $@

task7: $(BIN_DIR)/task7
$(BIN_DIR)/task7: $(OBJ_DIR)/task7.o $(MODULES)
	@mkdir -p $(BIN_DIR)
	$(CC_MPI) $(CFLAGS) $^ -o $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)