# Define directories for source and object files
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# Define the main variants (without extensions)
MAIN_VARIANTS = molecular_simulation benchmark_electrostatic benchmark_dispersion benchmark_bonds

# List all common source files (excluding main.cpp)
COMMON_SRC_FILES = $(filter-out $(SRC_DIR)/main.cpp, $(wildcard $(SRC_DIR)/**/*.cpp))

# List object files for common source files
COMMON_OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(COMMON_SRC_FILES))

# Compiler settings
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -pedantic -O2 -g -fopenmp
LDFLAGS = -fopenmp

# Targets for each variant
TARGETS = $(patsubst %, $(BIN_DIR)/%, $(MAIN_VARIANTS))

# Default target (compiles all variants)
all: $(TARGETS)

# Build each target separately
$(BIN_DIR)/%: $(COMMON_OBJ_FILES) $(OBJ_DIR)/%.o
	@mkdir -p $(BIN_DIR)
	$(CXX) $^ -o $@ $(LDFLAGS)

# Compile the common source files to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)  # Create the subdirectory in obj if it doesn't exist
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile the specific main.cpp variants
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up object and binary files
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Phony targets
.PHONY: all clean
