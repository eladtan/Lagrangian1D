CC := g++
DUMP_FOLDER := temp_$(TARGET)
PREXECUTE :=
ENDEMIC_SOURCES := $(shell find $(TARGET) -name '*.cpp')
ENDEMIC_OBJECTS := $(patsubst $(TARGET)/%.cpp,$(DUMP_FOLDER)/%.o,$(ENDEMIC_SOURCES))
ifeq ($(MODE),debug)
	OPTIMIZATION_FLAGS := -O0 -g -pg -std=c++11
	LFLAGS := -pg
else ifeq ($(MODE),intel)
	CC := icpc
	OPTIMIZATION_FLAGS := -O3 -ipo -xHost -fp-model precise -std=c++11
	LINT_FLAGS :=
	LFLAGS :=
	ARCHIVER_FUNC := xiar
else
	MODE := production
	OPTIMIZATION_FLAGS := -O3 -std=c++11
	LFLAGS :=
endif
CFLAGS := $(LINT_FLAGS) $(OPTIMIZATION_FLAGS)

$(DUMP_FOLDER)/normal_termination.res: $(DUMP_FOLDER)/test.exe

$(DUMP_FOLDER)/test.exe: $(ENDEMIC_OBJECTS) | ../library_$(MODE)/libLag.a
	echo "$(ENDEMIC_SOURCES)"
	$(CC) $^ $| -o $@ $(LFLAGS) -lhdf5 -lhdf5_cpp -L $(HDF5_LIB_PATH)

../library_$(MODE)/librich.a:
	$(MAKE) -C .. MODE=$(MODE)

$(ENDEMIC_OBJECTS): $(DUMP_FOLDER)/%.o: $(TARGET)/%.cpp
	mkdir -p $(DUMP_FOLDER)
	$(CC) $(CFLAGS) -c $< -o $@ -I ..
