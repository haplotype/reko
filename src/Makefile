# Makefile for compiling the sphere_vec_program with Eigen library

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -O2 -Wall -Wextra -pthread #-ggdb

# Include directories
EIGEN_DIR = /usr/local/include/eigen3/

# Target executable name
TARGET = reko

# Source files
SRC = knockoff.cpp thread_pool.cpp

# Object files
OBJ = $(SRC:.cpp=.o)

# Include directories
INCLUDES = -I$(EIGEN_DIR)

# Build rules
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<

clean:
	rm -f $(OBJ) $(TARGET)

