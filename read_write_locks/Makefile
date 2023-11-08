# This makes both a sequential and parallel program to display the Mandelbrot Set
# Use module load intel

CC = icc

CFLAGS_SERIAL = -g -O3 -xHost -fno-alias -std=c99 -I include

SRC_DIR 	= src
OBJ_DIR 	= obj
BIN_DIR 	= bin

SRC_FILES = $(shell find $(SRC_DIR) -name '*.c')
OBJ_FILES = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILES))

# For debugging
$(info SRC_FILES = $(SRC_FILES))
$(info OBJ_FILES = $(OBJ_FILES))

all: $(addprefix $(BIN_DIR)/, $(EXECUTABLES))

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS_SERIAL) -c $< -o $@

# Specific rule forlist 
list: $(BIN_DIR)/list
$(BIN_DIR)/list: $(OBJ_DIR)/list.o $(MODULES)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_SERIAL) $^ -o $@
