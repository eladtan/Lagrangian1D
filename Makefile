SOURCE_DIR := source
RAW_SOURCES := $(shell find $(SOURCE_DIR) -name '*.cpp')
SOURCES := $(RAW_SOURCES)
LIB_FILE = libLag.a
CC := g++
CCC := gcc
ARCHIVER_FUNC := ar
ifeq ($(MODE),debug)
	OPTIMIZATION_FLAGS := -O0 -g -pg -std=c++11
	LINT_FLAGS :=
else ifeq ($(MODE),intel)
	CC := icpc
	CCC := icc
	OPTIMIZATION_FLAGS := -O3 -ipo -xHost -fp-model precise -std=c++11
	LINT_FLAGS :=
	ARCHIVER_FUNC := xiar
else
	MODE = production
	OPTIMIZATION_FLAGS := -O3 -march=native -std=c++11 -fno-expensive-optimizations
endif
LIBRARY_FOLDER := library_$(MODE)
OBJECTS := $(patsubst $(SOURCE_DIR)/%.cpp,$(LIBRARY_FOLDER)/%.o,$(SOURCES))

$(LIBRARY_FOLDER)/$(LIB_FILE): $(OBJECTS)
	$(ARCHIVER_FUNC) cr $@ $^

$(OBJECTS): $(LIBRARY_FOLDER)/%.o: $(SOURCE_DIR)/%.cpp
	mkdir -p `dirname $@`
	$(CC) -c $(OPTIMIZATION_FLAGS) $(LINT_FLAGS) $< -o $@
	$(CC) -MM $(OPTIMIZATION_FLAGS) $(LINT_FLAGS) $< -o $(LIBRARY_FOLDER)/$*.d
	@sed 's,\(\w*\)\.o,$@,g' -i $(LIBRARY_FOLDER)/$*.d

-include $(OBJECTS:.o=.d)

clean:
	rm -rf ./$(LIBRARY_FOLDER)
