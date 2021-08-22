#
# Common Variables
#

CC      := g++
OBJ     := build
SRC      = src
SRCS     = $(wildcard $(SRC)/pathTracer/*.cpp)
SRCS    += $(wildcard $(SRC)/rayTracer/*.cpp)
SRCS    += $(wildcard $(SRC)/normalsDisplay/*.cpp)
SRCS    += $(wildcard $(SRC)/utils/*.cpp)
SRCS    += $(SRC)/Render3D.cpp
DEP     := deps
DEPS    := $(patsubst $(SRC)/%.cpp,$(DEP)/%.d,$(SRCS))
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEP)/$*.d
EXE     := Render3D

#
# Release Variables
#
OBJS     = $(patsubst $(SRC)/%.cpp,$(OBJ)/%.o,$(SRCS))
R_CFLAGS = -g -std=c++17 -O3 -ffast-math $(OPENMP) -ftree-vectorize -msse2 -mfpmath=sse -flto -march=native
R_LDLIBS = -lm -lpng $(OPENMP)

#
# Debug variables
#
DOBJ     = $(OBJ)_deb
DOBJS    = $(patsubst $(SRC)/%.cpp,$(DOBJ)/%.o,$(SRCS))
CFLAGS   = -g -std=c++17 -DDEBUG
LDLIBS   = -lm -lpng

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

$(DEPS):
	@mkdir -p "$(@D)"

# Build the executable
R$(EXE): $(OBJS)
	@echo "Building release executable " $(EXE)
	@$(CC) $(R_CFLAGS) $^ -o $(EXE) $(R_LDLIBS)

# Compile the individual CPP files
$(OBJ)/%.o: $(SRC)/%.cpp $(DEP)/%.d | $(DEP)
	@mkdir -p "$(@D)"
	@echo Compiling $<
	@$(CC) $(DEPFLAGS) $(R_CFLAGS) -c $< -o $@

# Create OBJ directory
$(OBJ):
	mkdir -p $@

#
# Debug rules
#

# No compiler optimizations enabled, use default flags.
debug: D$(EXE)

# Build the executable
D$(EXE): $(DOBJS)
	@echo "Building debug executable " $(EXE)
	@$(CC) $(CFLAGS) $^ -o $(EXE) $(LDLIBS)

# Compile the individual CPP files
$(DOBJ)/%.o: $(SRC)/%.cpp $(DEP)/%.d | $(DEP)
	@mkdir -p "$(@D)"
	@echo Compiling $<
	@$(CC) $(DEPFLAGS) $(CFLAGS) -c $< -o $@

# Create OBJ directory
$(DOBJ):
	mkdir -p $@

# Create Dependencies directory
$(DEP): 
	@mkdir -p $@

#
# Other rules
#

clean:
	rm -rf $(OBJ) $(DOBJ) $(EXE) $(DEP)

include $(wildcard $(DEPS))