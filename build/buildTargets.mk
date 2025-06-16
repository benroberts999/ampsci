# All makefile targets

################################################################################
# Executable targets
DEFAULTEXES = $(addprefix $(XD)/, ampsci)
ifneq ($(Build),release)
  DEFAULTEXES += $(XD)/tests
endif

ALLEXES = $(addprefix $(XD)/, ampsci tests)

CHECKS = GitInfo checkObj checkXdir

# Default make rule:
all: $(CHECKS) $(DEFAULTEXES)

# "Full" make rule: used for VSCode make/intellisence
full: $(CHECKS) $(ALLEXES)

################################################################################
# Automatically generate dependency files for each cpp file, + compile:
# Put built objects in subdirectories based on build mode and compiler

# Put built objects in sub-directory based on compiler and OMP setting
# OMP_BDIR is set in buildOptions.mk
SUBBD=$(Build)/$(notdir $(word 1, $(CXX)))/$(OMP_BDIR)

# Auto rule for all cpp files
$(BD)/$(SUBBD)/%.o: $(SD)/%.cpp
	@mkdir -p $(@D)
	$(COMP)

# Force version.o to build each time: ensure updated git+version info!
.PHONY: force
$(BD)/$(SUBBD)/version/version.o: $(SD)/version/version.cpp force
	@mkdir -p $(@D)
	$(COMP) $(GITFLAGS)

# include the dependency files
-include $(BD)/$(SUBBD)/*.d $(BD)/$(SUBBD)/*/*.d

################################################################################
# List all objects in sub-directories (i.e., that don't contain a main())
# Test files: all "test" files are named <filename>.tests.cpp
TEST_SRC_FILES := $(wildcard $(SD)/*/*.tests.cpp)
TEST_OBJS := $(subst $(SD),$(BD)/$(SUBBD),$(subst .cpp,.o,$(TEST_SRC_FILES)))

# Each non-test source file (except main())
SRC_FILES := $(filter-out $(TEST_SRC_FILES), $(wildcard $(SD)/*/*.cpp))
OBJS := $(subst $(SD),$(BD)/$(SUBBD),$(subst .cpp,.o,$(SRC_FILES)))

################################################################################
# Link + build all final programs

$(XD)/ampsci: $(BD)/$(SUBBD)/main.o $(OBJS) | $(XD)
	@echo $(BD)/$(SUBBD)/main.o $(OBJS) > $(BD)/objects.txt
	$(LINK)

$(XD)/tests: $(OBJS) $(TEST_OBJS) | $(XD)
	$(LINK)

################################################################################
# Directory creation rules

.PHONY: checkObj checkXdir

checkObj:
	@mkdir -p $(BD)

checkXdir:
	@mkdir -p $(XD)

################################################################################
# Git info and utility targets

NOW := $(shell date +%Y-%m-%d' '%H:%M' '%Z 2>/dev/null)
GITFLAGS = -D GITREVISION="$(shell git rev-parse --short HEAD 2>/dev/null)"
GITFLAGS += -D GITBRANCH="$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)"
GITFLAGS += -D GITMODIFIED="$(shell git status -s 2>/dev/null)"
GITFLAGS += -D CXXVERSION="$(shell $(CXX) --version 2>/dev/null | sed -n '1p' | sed s/'('/'['/ | sed s/')'/']'/)"
GITFLAGS += -D COMPTIME="$(NOW)"

.PHONY: GitInfo
GitInfo:
	@echo Git Branch: $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)
	@echo Git Revision: $(shell git rev-parse --short HEAD 2>/dev/null)
	@echo Modified: $(shell git status -s 2>/dev/null)

################################################################################
# Utility targets

.PHONY: all full clean remove_junk docs doxy do_the_chicken_dance clang_format pre_commit test includes checkObj checkXdir GitInfo force

clean:
	rm -f -v $(ALLEXES)
	rm -rf -v $(BD)/*.o $(BD)/*.d $(BD)/*.gc* $(BD)/*.txt $(BD)/*/
	rm -f -v *deleteme*

remove_junk:
	rm -f -v *deleteme*

# Builds the ampsci.pdf (physics) docs
docs:
	( cd ./doc/tex && make )
	cp ./doc/tex/ampsci.pdf ./doc/ampsci.pdf
	( cd ./doc/tex && make clean)

# Builds code documentation using doxygen, if available
doxy:
	doxygen ./doc/doxygen/Doxyfile
	( cd ./doc/latex && make )
	cp ./doc/latex/refman.pdf ./doc/documentation.pdf
	make docs
	cp ./doc/ampsci.pdf ./docs/ampsci.pdf 2>/dev/null || :
	( cd ./doc/latex && make clean)

do_the_chicken_dance:
	@echo 'Why would I do that?'

# Runs clang format on the entire project
CLANG_FORMAT_EXCLUDES := $(SD)/fmt $(SD)/catch2 $(SD)/json
clang_format:
	@echo "Running clang format (whole project)"
	find $(SD)/ \( $(foreach dir,$(CLANG_FORMAT_EXCLUDES),-path "$(dir)" -o ) -false \) -prune -o \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) -print | xargs -r clang-format -i -verbose

# Makes the includes
includes:
	./$(SD)/build_includes.sh
	@echo "Running clang format (on includes)"
	clang-format -i -verbose $(SD)/DiracOperator/Operators.hpp
	find $(SD)/ -iname 'include.hpp' | xargs clang-format -i -verbose

license-year:
	@year=$$(date +%Y); \
	sed -i.bak -E "s/(Copyright \(c\) )([0-9]{4})(–[0-9]{4})?/\1\2–$${year}/" LICENSE; \
	rm -f LICENSE.bak; \
	echo "Updated license year to $${year}"


pre_commit:
	make includes
	make clang_format
	make ampsci tests
	make test
	make remove_junk