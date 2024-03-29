# This makes both a sequential and parallel program to display the Mandelbrot Set
# Use module load intel

CC = icc

CFLAGS_OMP    = -g -O3 -xHost -fno-alias -std=c99 -qopenmp -I include 
CFLAGS_SERIAL = -g -O3 -xHost -fno-alias -std=c99 -I include

EXECUTABLES = mandseq mandseq-avx mandomp mandomp-avx mandomp-ts-avx mandomp-collapse mandomp-collapse-avx mandomp-tasks mandomp-tasks-row mandomp-tasks-row-shared

SRC_DIR 	= src
OBJ_DIR 	= obj
BIN_DIR 	= bin
INCLUDE_DIR = include

SRC_FILES = $(shell find $(SRC_DIR) -name '*.c')
OBJ_FILES = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILES))

DRAND  	  = $(OBJ_DIR)/drand.o
DRAND_TS  = $(OBJ_DIR)/drand-ts.o
TIMING 	  = $(OBJ_DIR)/timing.o
UTILITIES = $(OBJ_DIR)/utilities.o

MODULES    = $(DRAND) $(TIMING) $(UTILITIES)
MODULES_TS = $(DRAND_TS) $(TIMING) $(UTILITIES)

# For debugging
$(info SRC_FILES = $(SRC_FILES))
$(info OBJ_FILES = $(OBJ_FILES))

all: $(addprefix $(BIN_DIR)/, $(EXECUTABLES))

# General rule for object files (requiring OpenMP)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS_OMP) -c $< -o $@

# Specific rule for mandseq
$(OBJ_DIR)/mandseq.o: $(SRC_DIR)/mandseq.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS_SERIAL) -c $< -o $@	

$(OBJ_DIR)/mandseq-avx.o: $(SRC_DIR)/mandseq-avx.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS_SERIAL) -c $< -o $@	

# Specific rule for mandseq
mandseq: $(BIN_DIR)/mandseq
$(BIN_DIR)/mandseq: $(OBJ_DIR)/mandseq.o $(MODULES)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_SERIAL) $^ -o $@

mandseq-avx: $(BIN_DIR)/mandseq-avx
$(BIN_DIR)/mandseq-avx: $(OBJ_DIR)/mandseq-avx.o $(MODULES)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_SERIAL) $^ -o $@

# Specific rule for mandomp
mandomp: $(BIN_DIR)/mandomp
$(BIN_DIR)/mandomp: $(OBJ_DIR)/mandomp.o $(MODULES)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_OMP) $^ -o $@

mandomp-avx: $(BIN_DIR)/mandomp-avx
$(BIN_DIR)/mandomp-avx: $(OBJ_DIR)/mandomp-avx.o $(MODULES)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_OMP) $^ -o $@

# Specific rule for mandomp-ts
mandomp-ts: $(BIN_DIR)/mandomp-ts
$(BIN_DIR)/mandomp-ts: $(OBJ_DIR)/mandomp-ts.o $(MODULES_TS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_OMP) $^ -o $@

mandomp-ts-avx: $(BIN_DIR)/mandomp-ts-avx
$(BIN_DIR)/mandomp-ts-avx: $(OBJ_DIR)/mandomp-ts-avx.o $(MODULES_TS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_OMP) $^ -o $@

mandomp-collapse: $(BIN_DIR)/mandomp-collapse
$(BIN_DIR)/mandomp-collapse: $(OBJ_DIR)/mandomp-collapse.o $(MODULES_TS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_OMP) $^ -o $@

mandomp-collapse-avx: $(BIN_DIR)/mandomp-collapse-avx
$(BIN_DIR)/mandomp-collapse-avx: $(OBJ_DIR)/mandomp-collapse-avx.o $(MODULES_TS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_OMP) $^ -o $@

mandomp-tasks: $(BIN_DIR)/mandomp-tasks
$(BIN_DIR)/mandomp-tasks: $(OBJ_DIR)/mandomp-tasks.o $(MODULES_TS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_OMP) $^ -o $@

mandomp-tasks-columns: $(BIN_DIR)/mandomp-tasks-columns
$(BIN_DIR)/mandomp-tasks-columns: $(OBJ_DIR)/mandomp-tasks-columns.o $(MODULES_TS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_OMP) $^ -o $@

mandomp-tasks-columns-shared: $(BIN_DIR)/mandomp-tasks-columns-shared
$(BIN_DIR)/mandomp-tasks-columns-shared: $(OBJ_DIR)/mandomp-tasks-columns-shared.o $(MODULES_TS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_OMP) $^ -o $@

mandomp-ts-avx-parallel: $(BIN_DIR)/mandomp-ts-avx-parallel
$(BIN_DIR)/mandomp-ts-avx-parallel: $(OBJ_DIR)/mandomp-ts-avx-parallel.o $(MODULES_TS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS_OMP) $^ -o $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Specific clean targets for each executable
clean-mandseq: 
	rm -f $(BIN_DIR)/mandseq $(OBJ_DIR)/mandseq.o

# Specific clean targets for each executable
clean-mandseq-avx: 
	rm -f $(BIN_DIR)/mandseq-avx $(OBJ_DIR)/mandseq-avx.o

clean-mandomp: 
	rm -f $(BIN_DIR)/mandomp $(OBJ_DIR)/mandomp.o

clean-mandomp-avx: 
	rm -f $(BIN_DIR)/mandomp-avx $(OBJ_DIR)/mandomp-avx.o

clean-mandomp-ts:
	rm -f $(BIN_DIR)/mandomp-ts $(OBJ_DIR)/mandomp-ts.o

clean-mandomp-ts-avx:
	rm -f $(BIN_DIR)/mandomp-ts-avx $(OBJ_DIR)/mandomp-ts-avx.o

clean-mandomp-collapse:
	rm -f $(BIN_DIR)/mandomp-collapse $(OBJ_DIR)/mandomp-collapse.o

clean-mandomp-collapse-avx:
	rm -f $(BIN_DIR)/mandomp-collapse-avx $(OBJ_DIR)/mandomp-collapse-avx.o

clean-mandomp-tasks:
	rm -f $(BIN_DIR)/mandomp-tasks $(OBJ_DIR)/mandomp-tasks.o

clean-mandomp-tasks-columns:
	rm -f $(BIN_DIR)/mandomp-tasks-columns $(OBJ_DIR)/mandomp-tasks-columns.o

clean-mandomp-tasks-columns-shared:
	rm -f $(BIN_DIR)/mandomp-tasks-columns-shared $(OBJ_DIR)/mandomp-tasks-columns-shared.o

clean-mandomp-ts-avx-parallel:
	rm -f $(BIN_DIR)/mandomp-ts-avx-parallel $(OBJ_DIR)/mandomp-ts-avx-parallel.o