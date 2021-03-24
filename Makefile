ROOT_INC=$(shell root-config --cflags)
ROOT_LIB=$(shell root-config --libs)

LINKFLAGS += $(ROOT_LIB) -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -lpthread -lm -ldl
DEBUG =-O3
CPPFLAGS += $(ROOT_INC) -I./include -Wall -fpic -pipe -O3

SRC_DIR := src
INC_DIR := include
OBJ_DIR := obj
BIN_DIR := bin

SRC := $(wildcard $(SRC_DIR)/*.cc)
INC := $(wildcard $(INC_DIR)/*.h)
OBJ := $(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
EXE := $(BIN_DIR)/prog

.PHONY: clean all
all: $(EXE)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(INC_DIR)/%.h | $(OBJ_DIR)
	g++ $(CPPFLAGS) -c -o $@ $<

$(EXE): $(OBJ) | $(BIN_DIR)
	g++ $(LINKFLAGS) $(DEBUG) $^ -o $@

clean:
	-rm $(BIN_DIR)/prog $(OBJ_DIR)/*.o
