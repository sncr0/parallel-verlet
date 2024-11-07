# Define directories for source and object files
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# List all source files, including main.cpp
SRC_FILES = $(wildcard $(SRC_DIR)/main.cpp) $(wildcard $(SRC_DIR)/system/*.cpp)

# List all object files
OBJ_FILES = $(SRC_FILES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Compiler settings
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -pedantic -O2
LDFLAGS = 

# Output binary
TARGET = $(BIN_DIR)/molecular_simulation

# Default target (compiles everything)
all: $(TARGET)

# Link the object files into the binary
$(TARGET): $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(OBJ_FILES) -o $(TARGET) $(LDFLAGS)

# Compile the source files to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up object and binary files
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Phony targets
.PHONY: all clean
