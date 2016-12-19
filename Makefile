CORE = sift4g
DEPS = vendor/swsharp

all: TARGETS=cpu
gpu: TARGETS=all
clean: TARGETS=clean

all: $(CORE) $(DEPS)
gpu: $(CORE) $(DEPS)

clean: $(CORE) $(DEPS)
	@echo [RM] cleaning
	@rm $(OBJ_DIR) $(EXC_DIR) -rf

$(DEPS):
	@echo \>\>\> $@ \<\<\<
	@$(MAKE) -s -C $@ $(TARGETS)

$(CORE): $(DEPS)
	@echo \>\>\> $@ \<\<\<
	@$(MAKE) -s -C $@ $(TARGETS)

.PHONY: $(CORE) $(DEPS)
