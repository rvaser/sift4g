CORE = sift4g
DEPS = vendor/swsharp

all: TARGETS=all
cpu: TARGETS=cpu
clean: TARGETS=clean

all: $(CORE) $(DEPS)
cpu: $(CORE) $(DEPS)

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
