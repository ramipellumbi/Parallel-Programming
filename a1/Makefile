# Assumes that you have already loaded the intel compiler:
#        module load intel
#
# Assumes that pi.c is in the part1/ directory
# All C code in part1/ is built
# Assumes that timing.c and timing.h are in the common/ directory

CC = icc

# Various compiler options

# CFLAGS = -g -O0 -fno-alias -std=c99 -I common
# CFLAGS = -g -O1 -fno-alias -std=c99
# CFLAGS = -g -O0 -no-vec -no-simd -fno-alias -std=c99 -fp-model precise
# CFLAGS = -O0 -no-vec -no-simd -fno-alias -std=c99
# CFLAGS = -O1 -no-vec -no-simd -fno-alias -std=c99
# CFLAGS = -O3 -no-vec -no-simd -fno-alias -std=c99
# CFLAGS = -g -O3 -xHost -fno-alias -std=c99 -qopt-report=5 -fsource-asm
# CFLAGS = -std=c99

CFLAGS = -g -O3 -xHost -fno-alias -std=c99 -I common

# Build directories
BUILD_DIR := build
PART1_BUILD_DIR := $(BUILD_DIR)/part1
PART2_BUILD_DIR := $(BUILD_DIR)/part2
COMMON_BUILD_DIR := $(BUILD_DIR)/common

# Source code
COMMON_SRC := common/timing.c
PART1_SRC := $(wildcard part1/*.c)
PART2_SRC := $(wildcard part2/*.c)

# Object files
COMMON_OBJ := $(COMMON_BUILD_DIR)/timing.o
PART1_OBJ := $(PART1_SRC:.c=.o)
PART2_OBJ := $(PART2_SRC:.c=.o)

all: pi triad

$(COMMON_OBJ): $(COMMON_SRC)
	@mkdir -p $(COMMON_BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(PART1_BUILD_DIR)/%.o: part1/%.c
	@mkdir -p $(PART1_BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(PART2_BUILD_DIR)/%.o: part2/%.c
	@mkdir -p $(PART2_BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Build code for Exercise 1
pi: $(PART1_BUILD_DIR)/pi.o $(PART1_BUILD_DIR)/integrate.o $(COMMON_BUILD_DIR)/timing.o
	$(CC) -o $@ $(CFLAGS) $^

# Build code for Exercise 2
triad: $(PART2_BUILD_DIR)/triad.o $(PART2_BUILD_DIR)/dummy.o $(COMMON_BUILD_DIR)/timing.o
	$(CC) -o $@ $(CFLAGS) $^

clean:
	rm -rf $(BUILD_DIR) pi triad
