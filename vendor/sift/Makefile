CORE = sift
DEPS = blimps

all: TARGETS=all
cpu: TARGETS=all
clean: TARGETS=clean

all: $(CORE) $(DEPS)
cpu: all

clean: $(CORE) $(DEPS)
	@echo [RM] cleaning
	@rm $(OBJ_DIR) $(EXC_DIR) -rf

$(DEPS):
	@echo [DEP] $@
	@$(MAKE) -s -C $@ $(TARGETS)

$(CORE): $(DEPS)
	@echo [CORE] $@
	@$(MAKE) -s -C $@ $(TARGETS)

.PHONY: $(CORE) $(DEPS)
