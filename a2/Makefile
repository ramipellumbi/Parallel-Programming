# This makes both a sequential and parallel program to display the Mandelbrot Set
# Use module load intel

CC = icc

CFLAGS_OMP 	  = -g -O3 -xHost -fno-alias -std=c99 -qopenmp -I include 
CFLAGS_SERIAL = -g -O3 -xHost -fno-alias -std=c99 -I include 

EXECUTABLES = mandseq mandomp mandomp-ts mandomp-collapse mandomp-tasks mandomp-tasks-row mandomp-tasks-row-shared

SRC_DIR 	= src
OBJ_DIR 	= obj
BIN_DIR 	= bin
INCLUDE_DIR = include

MAIN_TARGET = $(BIN_DIR)/main

SRC_FILES = $(shell find $(SRC_DIR) -name '*.c')
OBJ_FILES = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILES))

all: $(MAIN_TARGET)

$(MAIN_TARGET): $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_SERIAL) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS_SERIAL) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean
