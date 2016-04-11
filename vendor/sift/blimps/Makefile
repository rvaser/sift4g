CC = gcc
LD = gcc

NAME = blimps

OBJ_DIR = obj
SRC_DIR = src
INC_DIR = include
LIB_DIR = lib

I_CMD = $(addprefix -I, $(SRC_DIR))
L_CMD = $(addprefix -L, )

CC_FLAGS = $(I_CMD) --pedantic -std=c99 -O2 -Wall -Wno-unused-result -Wno-strict-overflow -Wno-aggressive-loop-optimizations
LD_FLAGS = $(I_CMD) $(L_CMD)

API = $(addprefix $(SRC_DIR)/, global.h errors.h memory.h residues.h blocks.h \
    matrix.h sequences.h aabet.h ntbet.h convert.h frequency.h files.h \
	protomat.h blocksprogs.h license.h strutil.h blastapp.h alphabet.h gcode.h \
	output.h pattern.h)

SRC = $(shell find $(SRC_DIR) -type f -regex ".*\.c" -not -name "blimps\.c")
OBJ = $(subst $(SRC_DIR), $(OBJ_DIR), $(addsuffix .o, $(basename $(SRC))))
DEP = $(OBJ:.o=.d)
INC = $(subst $(SRC_DIR), $(INC_DIR)/$(NAME), $(API))
LIB = $(LIB_DIR)/lib$(NAME).a

all: lib include
install: all

include: $(INC)

lib: $(LIB)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@echo [CC] $<
	@mkdir -p $(dir $@)
	@$(CC) $< -c -o $@ -MMD $(CC_FLAGS)

$(INC_DIR)/$(NAME)/%.h: $(SRC_DIR)/%.h
	@echo [CP] $@
	@mkdir -p $(dir $@)
	@cp $< $@

$(LIB): $(OBJ)
	@echo [AR] $@
	@mkdir -p $(dir $@)
	@ar rcs $(LIB) $(OBJ)

clean:
	@echo [RM] cleaning
	@rm $(LIB_DIR) $(INC_DIR) $(OBJ_DIR) -rf

-include $(DEP)
