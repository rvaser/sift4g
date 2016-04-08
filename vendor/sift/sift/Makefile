CC = gcc
LD = gcc

NAME =  sift

OBJ_DIR = obj
SRC_DIR = src
DOC_DIR = doc
INC_DIR = ../include/$(NAME)
LIB_DIR = ../lib
EXC_DIR = ../bin

BLIMPS_DIR = ../blimps

I_CMD = $(addprefix -I, $(SRC_DIR) $(BLIMPS_DIR)/include)
L_CMD = $(addprefix -L, $(BLIMPS_DIR)/lib)

DEP_LIBS = $(BLIMPS_DIR)/lib/libblimps.a

CC_FLAGS = $(I_CMD) -O2 -w
LD_FLAGS = $(I_CMD) $(L_CMD) -lblimps -lm

#API = $(addprefix $(SRC_DIR)/, )

SRC = $(shell find $(SRC_DIR) -type f -regex ".*\.c")
OBJ = $(subst $(SRC_DIR), $(OBJ_DIR), $(addsuffix .o, $(basename $(SRC))))
DEP = $(OBJ:.o=.d)
INC = $(subst $(SRC_DIR), $(INC_DIR), $(API))
LIB = $(LIB_DIR)/lib$(NAME).a
EXC = $(NAME)
BIN = $(EXC_DIR)/$(EXC)

all: $(EXC)

#install: lib include

#include: $(INC)

#lib: $(LIB)

$(EXC): $(OBJ) $(DEP_LIBS)
	@echo [LD] $@
	@mkdir -p $(dir $@)
	@$(LD) -o $@ $^ $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@echo [CC] $<
	@mkdir -p $(dir $@)
	@$(CC) $< -c -o $@ -MMD $(CC_FLAGS)

$(INC_DIR)/%.h: $(SRC_DIR)/%.h
	@echo [CP] $@
	@mkdir -p $(dir $@)
	@cp $< $@

$(LIB): $(OBJ)
	@echo [AR] $@
	@mkdir -p $(dir $@)
	@ar rcs $(LIB) $(OBJ)

$(BIN): $(EXC)
	@echo [CP] $@
	@mkdir -p $(dir $@)
	@cp $< $@

clean:
	@echo [RM] cleaning
	@rm $(OBJ_DIR) $(EXC) -rf

-include $(DEP)
