CC       = g++
OBJ      = build
SRC      = src
SRCS     = $(wildcard $(SRC)/pathTracer/*.cpp)
SRCS    += $(wildcard $(SRC)/rayTracer/*.cpp)
SRCS    += $(wildcard $(SRC)/utils/*.cpp)
SRCS    += $(SRC)/Render3D.cpp
OBJS     = $(patsubst $(SRC)/%.cpp,$(OBJ)/%.o,$(SRCS))
EXE      = Render3D
CFLAGS   = -g -I$(INCLUDE) -std=c++17
LDLIBS   = -lm

.PHONY: all run clean release debug

all: release

# Add correct OMP flags for Linux/MacOS
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	OPENMP += -fopenmp
endif
ifeq ($(UNAME_S),Darwin)
	OPENMP += -Xpreprocessor -fopenmp
	LDLIBS += -lomp
endif

# Build the executable
$(EXE): $(OBJS) | $(BIN)
	@echo "Building executable " $@
	@$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

# Compile the individual CPP files
$(OBJ)/%.o: $(SRC)/%.cpp | $(OBJ)
	@mkdir -p "$(@D)"
	@echo Compiling $<
	@$(CC) $(CFLAGS) -c $< -o $@

# Create OBJ directory
$(OBJ):
	mkdir -p $@

# No compiler optimizations enabled, use default flags.
debug: $(EXE)

# Optimize in release mode
release: CFLAGS += -O3 -ffast-math $(OPENMP) -ftree-vectorize -msse2 -mfpmath=sse -flto=full -march=native
release: LDLIBS += $(OPENMP)
release: $(EXE)

clean:
	rm -rf $(OBJ) $(EXE)

