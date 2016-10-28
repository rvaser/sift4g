CP = g++
LD = nvcc

NAME = sift4g

SRC_DIR = src
OBJ_DIR = obj
EXC_DIR = ../bin
VND_DIR = ../vendor

I_CMD = $(addprefix -I, $(SRC_DIR) $(VND_DIR)/swsharp/include $(VND_DIR)/sift/include)
L_CMD = $(addprefix -L, $(VND_DIR)/swsharp/lib $(VND_DIR)/sift/lib)

DEP_LIBS = $(VND_DIR)/swsharp/lib/libswsharp.a $(VND_DIR)/sift/lib/libsift.a

CP_FLAGS = $(I_CMD) -std=c++11 -O3 -Wall -Wno-write-strings
LD_FLAGS = $(I_CMD) $(L_CMD) -lswsharp -lsift -lblimps -lpthread -lm

SRC = $(shell find $(SRC_DIR) -type f -regex ".*\.cpp")
OBJ = $(subst $(SRC_DIR), $(OBJ_DIR), $(addsuffix .o, $(basename $(SRC))))
DEP = $(OBJ:.o=.d)
BIN = $(EXC_DIR)/$(NAME)

cpu: LD = $(CP)

all: $(BIN)
cpu: all

clean:
	@echo [RM] cleaning
	@rm $(EXC_DIR) $(OBJ_DIR) -rf

$(BIN): $(OBJ) $(DEP_LIBS)
	@echo [LD] $@
	@mkdir -p $(dir $@)
	@$(LD) $(OBJ) -o $@ $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo [CP] $<
	@mkdir -p $(dir $@)
	@$(CP) $< -c -o $@ -MMD $(CP_FLAGS)

-include $(DEP)
