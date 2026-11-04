CXX = g++

CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -Iinclude

LIB_NAME = libKaryoKE.a
TARGET = karyoke

SRC_DIR = source
OBJ_DIR = objects
INC_DIR = include
MAIN_SRC = $(SRC_DIR)/main.cpp
MAIN_OBJ = $(OBJ_DIR)/main.o

# Recursively find all files matching a pattern under a directory,
# at any depth
rwildcard = $(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))


ALL_SOURCES = $(call rwildcard,$(SRC_DIR)/,*.cpp)


# Library sources = everything except main.cpp, main.cpp is the CLI
# entry point, not part of the reusable library.
LIB_SOURCES = $(filter-out $(MAIN_SRC),$(ALL_SOURCES))

LIB_OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(LIB_SOURCES))

HEADERS = $(call rwildcard,$(INC_DIR)/,*.h)


.PHONY: all lib clean

all: $(TARGET)

# --------------------------------------------------
# Static library - everything except main.cpp
# --------------------------------------------------
lib: $(LIB_NAME)

$(LIB_NAME): $(LIB_OBJECTS)
	ar rcs $(LIB_NAME) $(LIB_OBJECTS)

# --------------------------------------------------
# CLI executable links against the static library,
# exactly as an external user's program would.
# --------------------------------------------------
$(TARGET): $(MAIN_OBJ) $(LIB_NAME)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(MAIN_OBJ) $(LIB_NAME) -lcairo


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -rf $(OBJ_DIR)
	rm -f $(TARGET) $(LIB_NAME)
