#
# Common Variables
#

CC      := g++
OBJ     := build
SRC      = src
SRCS     = $(wildcard $(SRC)/pathTracer/*.cpp)
SRCS    += $(wildcard $(SRC)/rayTracer/*.cpp)
SRCS    += $(wildcard $(SRC)/utils/*.cpp)
SRCS    += $(SRC)/Render3D.cpp
EXE     := Render3D

#
# Release Variables
#
OBJS     = $(patsubst $(SRC)/%.cpp,$(OBJ)/%.o,$(SRCS))
R_CFLAGS = -g -I$(INCLUDE) -std=c++17 -O3 -ffast-math $(OPENMP) -ftree-vectorize -msse2 -mfpmath=sse -flto -march=native
R_LDLIBS = -lm $(OPENMP)

#
# Debug variables
#
DOBJ     = $(OBJ)_deb
DOBJS    = $(patsubst $(SRC)/%.cpp,$(DOBJ)/%.o,$(SRCS))
CFLAGS   = -g -I$(INCLUDE) -std=c++17
LDLIBS   = -lm

.PHONY: all run clean release debug

all: release

# Add correct OMP flags for Windows/Linux/MacOS
ifeq ($(OS),Windows_NT)
	OBJ := $(OBJ)_Win
else
    UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		OBJ := $(OBJ)_Lin
		OPENMP += -fopenmp
	endif
	ifeq ($(UNAME_S),Darwin)
		OBJ := $(OBJ)_Mac
		OPENMP += -Xpreprocessor -fopenmp
		R_LDLIBS += -lomp
	endif
endif

#
# Release rules
#

# Optimize in release mode
release: R$(EXE)

# Build the executable
R$(EXE): $(OBJS) | $(BIN)
	@echo "Building release executable " $(EXE)
	@$(CC) $(R_CFLAGS) $^ -o $(EXE) $(R_LDLIBS)

# Compile the individual CPP files
$(OBJ)/%.o: $(SRC)/%.cpp | $(OBJ)
	@mkdir -p "$(@D)"
	@echo Compiling $<
	@$(CC) $(R_CFLAGS) -c $< -o $@

# Create OBJ directory
$(OBJ):
	mkdir -p $@

#
# Debug rules
#

# No compiler optimizations enabled, use default flags.
debug: D$(EXE)

# Build the executable
D$(EXE): $(DOBJS) | $(BIN)
	@echo "Building debug executable " $(EXE)
	@$(CC) $(CFLAGS) $^ -o $(EXE) $(LDLIBS)

# Compile the individual CPP files
$(DOBJ)/%.o: $(SRC)/%.cpp | $(DOBJ)
	@mkdir -p "$(@D)"
	@echo Compiling $<
	@$(CC) $(CFLAGS) -c $< -o $@

# Create OBJ directory
$(DOBJ):
	mkdir -p $@

#
# Other rules
#

clean:
	rm -rf $(OBJ) $(DOBJ) $(EXE)
