# BUILD TARGETS:

################################################################################
$(info )
$(info ============================================================================)
$(info AMPSCI Build Configuration:)
$(info   Compiled      : $(NOW))
$(info   Mode          : $(MODE))
$(info   Compiler      : $(CXXPATH) [$(CXXVERSION)])
$(info   CXXSTD        : $(CXXSTD))
$(info   Optimisation  : $(OPT))
$(info   OpenMP        : $(if $(OMPLIB),$(OMPLIB),disabled))
$(info   Libraries     : $(LDLIBS))
$(info   Build dir     : $(BUILD_DIR))
$(info   git revision  : $(GITBRANCH)/$(GITREVISION))
$(if $(strip $(GITMODIFIED)),$(info   Modified      : $(GITMODIFIED)))
$(info ============================================================================)
$(info )

################################################################################
## Default targets (depends on build mode)

DEFAULTEXES := ampsci
ifneq ($(MODE),release)
  DEFAULTEXES += tests
endif

ALLEXES := ampsci tests

.DEFAULT_GOAL := all

all:  $(DEFAULTEXES)
full: $(ALLEXES)

################################################################################
## Main targets

## Compile+link main ampsci
ampsci: $(BUILD_DIR)/main.o $(LIB_OBJS) | $(BUILD_DIR)
	$(LINK)

## Compile+link tests
tests: $(BUILD_DIR)/tests.o $(LIB_OBJS) $(TEST_OBJS) | $(BUILD_DIR)
	$(LINK)

## Force version.o to rebuild each time (_any_ changes should be reflected)
.PHONY: force
$(BUILD_DIR)/version/version.o: $(SRC)/version/version.cpp force | $(BUILD_DIR)
	@mkdir -p $(@D)
	$(COMPILE) $(GITFLAGS)

################################################################################
## Convenience targets:

.PHONY: all full clean remove_junk docs open clang_format pre_commit includes license-year force do_the_chicken_dance

clean:
	rm -fv $(ALLEXES)
	rm -rfv $(BUILD)/* $(BUILD)/*/

## Remove unit test outputs
remove_junk:
	rm -f -v *deleteme*

## Builds the ampsci.pdf docs using tex, and the html using doxygen, if available
docs:
	bash ./doc/build_docs.sh

do_the_chicken_dance:
	@echo 'Why would I do that?'

## clang-format targets
## External libraries (fmt, catch2, json) are excluded via their own .clang-format files
CLANG_FORMAT := clang-format

## Format the entire project (asks for confirmation)
clang_format_all:
	@./$(SRC)/tools/clang_format_all.sh $(CLANG_FORMAT) $(SRC)

## Format only files changed since last commit (git diff HEAD)
clang_format:
	@./$(SRC)/tools/clang_format.sh $(CLANG_FORMAT)

## Makes the include.hpp files, runs license-year and clang_format
includes:
	./$(SRC)/tools/build_includes.sh
	./$(SRC)/tools/license_year.sh
	./$(SRC)/tools/clang_format.sh $(CLANG_FORMAT)

## Updates licence year
license-year:
	@./$(SRC)/tools/license_year.sh

## All things should be done before major commits
pre_commit: includes clang_format license-year ampsci tests
	./test [unit]
	rm -f -v *deleteme*
