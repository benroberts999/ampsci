## BUILD RULES:

## Use subdirectories for Build, to allow faster compiler switches
ifeq ($(OMPLIB),)
  OMPDIR := noomp
else
  OMPDIR := omp
endif

BUILD_DIR := $(BUILD)/$(MODE)/$(CXXNAME)/$(OPT)/$(OMPDIR)

################################################################################
## Find source files and objects

## Discover all .cpp sources recursively
SRCS := $(patsubst $(SRC)/%,%,$(shell find $(SRC) -type f -name '*.cpp' -print))

## And the subset of which are "tests":
TEST_SRCS := $(patsubst $(SRC)/%,%,$(shell find $(SRC) -type f -name '*.tests.cpp' -print))

## Main sources are in src/ (not in subdirectory)
MAIN_SRCS := $(notdir $(wildcard $(SRC)/*.cpp))

## Library sources are any of the others (not mains, not tests)
LIB_SRCS  := $(filter-out $(MAIN_SRCS) $(TEST_SRCS),$(SRCS))

## Map source -> object (preserve relative path under $(SRC))
MAIN_OBJS := $(MAIN_SRCS:%.cpp=$(BUILD_DIR)/%.o)
LIB_OBJS  := $(LIB_SRCS:%.cpp=$(BUILD_DIR)/%.o)
TEST_OBJS := $(TEST_SRCS:%.cpp=$(BUILD_DIR)/%.o)

################################################################################
## Dependencies

## Dependency files (auto generated)
DEPS      := $(MAIN_OBJS:.o=.d) $(LIB_OBJS:.o=.d) $(TEST_OBJS:.o=.d)

## Dependency includes
-include $(DEPS)

################################################################################
## Rules for compilation:

## Ensure build directory exists
$(BUILD_DIR):
	@mkdir -p $@

## Compile all objects
$(BUILD_DIR)/%.o: $(SRC)/%.cpp | $(BUILD_DIR)
	@mkdir -p $(@D)
	$(COMPILE)

################################################################################
## External modules (set EXTERNAL_MODULES to a space-separated list of .cpp files)

## Function: maps a source path to its object file under ExternalModules/
## Uses basename only - all external module filenames must be unique.
ext_obj = $(BUILD_DIR)/ExternalModules/$(notdir $(1:.cpp=)).o

## Expand globs and resolve wildcards from EXTERNAL_MODULES (supports e.g. ext/*.cpp)
EXT_SRCS := $(filter %.cpp,$(wildcard $(EXTERNAL_MODULES)))

ifneq ($(strip $(EXT_SRCS)),)
  ## Compute object paths for all external modules
  EXT_OBJS := $(foreach src,$(EXT_SRCS),$(call ext_obj,$(src)))
  ## Include compiler-generated dependency files (.d) so header changes trigger recompilation
  -include $(EXT_OBJS:.o=.d)

  ## Rule template: generates one compile rule per external module.
  ## Called via $(eval $(call ext_rule,<src>)) below; $$(COMPILE) and $$(@D)
  ## use $$ so they survive the call expansion and are evaluated at recipe run time.
  define ext_rule
$(call ext_obj,$(1)): $(abspath $(1)) | $(BUILD_DIR)
	@mkdir -p $$(@D)
	$$(COMPILE)
  endef
  ## Instantiate a compile rule for each external source file
  $(foreach src,$(EXT_SRCS),$(eval $(call ext_rule,$(src))))
endif

