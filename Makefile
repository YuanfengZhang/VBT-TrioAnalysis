CC := g++
# CC := clang --analyze # and comment out the linker last line for sanity

BUILDDIR := build
TARGET := vbt


INCCORE := Core/include
INCDUO := DuoComparison/include
INCTRIO := MendelianViolation/include
INCVCFIO := VcfIO/include
INCUTIL := Utils
INCBASE := Base

CFLAGS := -std=c++11 -Wall -O2 -g
HTSLIB_PREFIX ?=
HTSINC := $(if $(HTSLIB_PREFIX),-I $(HTSLIB_PREFIX)/include -I $(HTSLIB_PREFIX)/include/htslib,-I htslib)
comma := ,
HTSLIB := $(if $(HTSLIB_PREFIX),-L $(HTSLIB_PREFIX)/lib -Wl$(comma)-rpath$(comma)$(HTSLIB_PREFIX)/lib,) -lhts
LIB := -lz -pthread $(HTSLIB)
INC := $(HTSINC) -I $(INCCORE) -I $(INCDUO) -I $(INCTRIO) -I $(INCVCFIO) -I $(INCUTIL) -I $(INCBASE) -I $(shell pwd)

SRCCORE := Core/src
SRCDUO := DuoComparison/src
SRCTRIO := MendelianViolation/src
SRCVCFIO := VcfIO/src
SRCUTIL := Utils
SRCBASE := Base
 
SOURCESCORE := $(shell find $(SRCCORE) -type f -name '*.cpp')
SOURCESDUO := $(shell find $(SRCDUO) -type f -name '*.cpp')
SOURCESTRIO := $(shell find $(SRCTRIO) -type f -name '*.cpp')
SOURCESVCFIO := $(shell find $(SRCVCFIO) -type f -name '*.cpp')


OBJECTSCORE := $(subst $(SRCCORE), $(BUILDDIR), $(SOURCESCORE:.cpp=.o))
OBJECTSDUO := $(subst $(SRCDUO), $(BUILDDIR), $(SOURCESDUO:.cpp=.o))
OBJECTSTRIO := $(subst $(SRCTRIO), $(BUILDDIR), $(SOURCESTRIO:.cpp=.o))
OBJECTSVCFIO := $(subst $(SRCVCFIO), $(BUILDDIR), $(SOURCESVCFIO:.cpp=.o))
OBJECTSUTIL := $(BUILDDIR)/CUtils.o
OBJECTSBASE := $(BUILDDIR)/CBaseVariantProvider.o

OBJECTS := $(OBJECTSCORE) $(OBJECTSDUO) $(OBJECTSTRIO) $(OBJECTSVCFIO) $(OBJECTSUTIL) $(OBJECTSBASE) $(BUILDDIR)/main.o

all: $(TARGET)
	@echo "SUCCESSFULLY COMPILED!!"
	

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " OBJECTS CORE: $(OBJECTSCORE)"
	@echo " OBJECTS DUO : $(OBJECTSDUO)"
	@echo " OBJECTS TRIO: $(OBJECTSTRIO)"
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB) $(INC)

$(BUILDDIR)/%.o: $(SRCCORE)/%.cpp Constants.h
	@mkdir -p $(BUILDDIR)
	@echo " CORE: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCDUO)/%.cpp Constants.h
	@mkdir -p $(BUILDDIR)
	@echo " DUO: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCTRIO)/%.cpp Constants.h
	@mkdir -p $(BUILDDIR)
	@echo " TRIO: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCVCFIO)/%.cpp Constants.h
	@mkdir -p $(BUILDDIR)
	@echo " VCFIO: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCBASE)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " BASE: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/%.o: $(SRCUTIL)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " UTILS: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@

$(BUILDDIR)/main.o: main.cpp
	@mkdir -p $(BUILDDIR)
	@echo " MAIN: $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c $< -o $@
 
clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)


.PHONY: clean

$(BUILDDIR)/mixed_replay_test: tests/mixed_replay.cpp $(filter-out $(BUILDDIR)/main.o,$(OBJECTS))
	$(CC) $(CFLAGS) $(INC) $^ -o $@ $(LIB)

test-core: $(BUILDDIR)/mixed_replay_test
	./$(BUILDDIR)/mixed_replay_test

.PHONY: test-core
