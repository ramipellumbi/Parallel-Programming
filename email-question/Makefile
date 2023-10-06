# This makes both a sequential and parallel program to display the Mandelbrot Set
# Use module load intel

CC = icc

CFLAGS_OMP    = -g -O3 -xHost -fno-alias -std=c99 -qopenmp -I include 
CFLAGS_SERIAL = -g -O3 -xHost -fno-alias -std=c99 -I include

EXECUTABLES = mandseq-avx 

SRC_DIR 	= src
OBJ_DIR 	= obj
BIN_DIR 	= bin
INCLUDE_DIR = include

SRC_FILES = $(shell find $(SRC_DIR) -name '*.c')
OBJ_FILES = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILES))

DRAND  	  = $(OBJ_DIR)/drand.o
TIMING 	  = $(OBJ_DIR)/timing.o

MODULES = $(DRAND) $(TIMING)

# For debugging
$(info SRC_FILES = $(SRC_FILES))
$(info OBJ_FILES = $(OBJ_FILES))

all: $(addprefix $(BIN_DIR)/, $(EXECUTABLES))

# General rule for object files (not requiring OpenMP)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS_SERIAL) -c $< -o $@

mandseq-avx: $(BIN_DIR)/mandseq-avx
$(BIN_DIR)/mandseq-avx: $(OBJ_DIR)/mandseq-avx.o $(MODULES)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_SERIAL) $^ -o $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)
