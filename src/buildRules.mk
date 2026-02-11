# BUILD RULES:

# Use subdirectories for Build, to allow faster compiler switches
ifeq ($(strip $(OMPLIB)),)
  OMPDIR := noomp
else
  OMPDIR := omp
endif

BUILD_DIR := $(BUILD)/$(MODE)/$(CXXNAME)/$(OPT)/$(OMPDIR)

################################################################################
# Rules for compiliation:

# Discover all .cpp sources recursively
SRCS := $(patsubst $(SRC)/%,%,$(shell find $(SRC) -type f -name '*.cpp' -print))

# And the subset of which are "tests":
TEST_SRCS := $(patsubst $(SRC)/%,%,$(shell find $(SRC) -type f -name '*.tests.cpp' -print))

# Main sources are in src/ (not in subdirectory)
MAIN_SRCS := $(notdir $(wildcard $(SRC)/*.cpp))

# Library sources are any of the others (not mains, not tests)
LIB_SRCS  := $(filter-out $(MAIN_SRCS) $(TEST_SRCS),$(SRCS))

# ------------------------------------------------------------
# Map source -> object (preserve relative path under $(SRC))

MAIN_OBJS := $(MAIN_SRCS:%.cpp=$(BUILD_DIR)/%.o)
LIB_OBJS  := $(LIB_SRCS:%.cpp=$(BUILD_DIR)/%.o)
TEST_OBJS := $(TEST_SRCS:%.cpp=$(BUILD_DIR)/%.o)

DEPS := $(MAIN_OBJS:.o=.d) $(LIB_OBJS:.o=.d) $(TEST_OBJS:.o=.d)


$(BUILD_DIR):
	@mkdir -p $@

# Compile all objects
$(BUILD_DIR)/%.o: $(SRC)/%.cpp | $(BUILD_DIR)
	@mkdir -p $(@D)
	$(COMPILE)

ampsci: $(BUILD_DIR)/main.o $(LIB_OBJS) | $(BUILD_DIR)
	$(LINK)

tests: $(BUILD_DIR)/tests.o $(LIB_OBJS) $(TEST_OBJS) | $(BUILD_DIR)
	$(LINK)

# Dependency includes
-include $(DEPS)

################################################################################
# Git info for version details

# These are passed as CMD-line input to the code, so code compiles with git info!
NOW := $(shell date +%Y-%m-%d' '%H:%M' '%Z 2>/dev/null)
GITFLAGS = -D GITREVISION="$(shell git rev-parse --short HEAD 2>/dev/null)"
GITFLAGS += -D GITBRANCH="$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)"
GITFLAGS += -D GITMODIFIED="$(shell git status -s 2>/dev/null)"
GITFLAGS += -D CXXVERSION="$(shell $(CXX) --version 2>/dev/null | sed -n '1p' | sed s/'('/'['/ | sed s/')'/']'/)"
GITFLAGS += -D COMPTIME="$(NOW)"

# Force version.o to rebuild each time
.PHONY: force
$(BUILD_DIR)/version/version.o: $(SRC)/version/version.cpp force | $(BUILD_DIR)
	@mkdir -p $(@D)
	$(COMPILE) $(GITFLAGS)
