# All makefile targets

################################################################################
#Allow exectuables to be placed in another directory:
DEFAULTEXES = $(addprefix $(XD)/, \
 ampsci \
)
# Compile 'tests' by default, unless in 'release' build
ifneq ($(Build),release)
  DEFAULTEXES += $(XD)/tests
endif

ALLEXES = $(addprefix $(XD)/, \
 ampsci tests \
)

CHECKS = GitInfo checkObj checkXdir

# Default make rule:
all: $(CHECKS) $(DEFAULTEXES)

# "Full" make rule: used for VSCode make/intellisence
full: $(CHECKS) $(ALLEXES)

################################################################################
# Automatically generate dependency files for each cpp file, + compile:
# Put built objects in subdirectories based on build mode and compiler

OMP_BDIR=omp
ifneq ($(UseOpenMP),true)
ifneq ($(UseOpenMP),yes)
  OMP_BDIR=
endif
endif

# Put built objects in sub-directory based on compiler and OMP seeting
# This is mainly for testing: faster re-compiles when switching settings!
SUBBD=$(Build)/$(notdir $(word 1, $(CXX)))/$(OMP_BDIR)

# Auto rule for all cpp files
$(BD)/$(SUBBD)/%.o: $(SD)/%.cpp
	@mkdir -p $(@D)
	$(COMP)

# Force version.o to build each time: ensure updated git+version info!
# Not the best solution, but works reasonably well
.PHONY: force
$(BD)/$(SUBBD)/version/version.o: $(SD)/version/version.cpp force
	@mkdir -p $(@D)
	$(COMP) $(GITFLAGS)

# include the dependency files
-include $(BD)/$(SUBBD)/*.d $(BD)/$(SUBBD)/*/*.d

################################################################################
# List all objects in sub-directories (i.e., that don't conatain a main())
# Automatically create list of target objects based on string matching.
# All main() files are: src/*.cpp
# All non-main files are: src/sub-dir/*cpp
# e.g., for each src/sub_dir/file.cpp -> build/sub_dir/file.o

# Test files: all "test" files are names <filename>.tests.cpp
TEST_SRC_FILES := $(wildcard $(SD)/*/*.tests.cpp)
TEST_OBJS := $(subst $(SD),$(BD)/$(SUBBD),$(subst .cpp,.o,$(TEST_SRC_FILES)))

# Each non-test source file (except main())
SRC_FILES := $(filter-out $(TEST_SRC_FILES), $(wildcard $(SD)/*/*.cpp))
OBJS := $(subst $(SD),$(BD)/$(SUBBD),$(subst .cpp,.o,$(SRC_FILES)))

################################################################################
# Link + build all final programs

$(XD)/ampsci: $(BD)/$(SUBBD)/main.o $(OBJS)
	@echo $(BD)/$(SUBBD)/main.o $(OBJS) > $(BD)/objects.txt
	$(LINK)

$(XD)/tests: $(OBJS) $(TEST_OBJS)
	$(LINK)

################################################################################

# Uses command-line parsing to get info from git and compiler.
# Gets: commit hash, branch name, which files have been modified since 
# last commit, the c++ compiler version, and date/time of compilation.
# It's set up so that it won't fail if git isn't installed (will be blank).
# It then passes this info via command-line arguments (-D) to compiler.
NOW:=$(shell date +%Y-%m-%d' '%H:%M' '%Z 2>/dev/null)
GITFLAGS=-D GITREVISION="$(shell git rev-parse --short HEAD 2>/dev/null)"
GITFLAGS+=-D GITBRANCH="$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)"
GITFLAGS+=-D GITMODIFIED="$(shell git status -s 2>/dev/null)"
GITFLAGS+=-D CXXVERSION="$(shell $(CXX) --version 2>/dev/null | sed -n '1p' | sed s/'('/'['/ | sed s/')'/']'/)"
GITFLAGS+=-D COMPTIME="$(NOW)"

################################################################################
################################################################################

# Just prints some git information to screen - has no impact
GitInfo:
	@echo Git Branch: $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)
	@echo Git Revision: $(shell git rev-parse --short HEAD 2>/dev/null)
	@echo Modified: $(shell git status -s 2>/dev/null)

# Check to see if 'obj' sub-directory exists. Doesn't make it.
checkObj:
	@if [ ! -d $(BD) ]; then \
	  echo '\n ERROR: Directory: '$(BD)' doesnt exist - please create it!\n'; \
	  false; \
	fi

# Check to see if 'excecutable' sub-directory exists. Doesn't make it.
checkXdir:
	@if [ ! -d $(XD) ]; then \
		echo '\n ERROR: Directory: '$(XD)' doesnt exist - please create it!\n'; \
		false; \
	fi

# Makes the includes
includes:
	./$(SD)/build_includes.sh
	@echo "Running clang format (on includes)"
	clang-format -i -verbose $(SD)/DiracOperator/Operators.hpp
	find $(SD)/ -iname 'include.hpp' | xargs clang-format -i -verbose

.PHONY: clean docs doxy do_the_chicken_dance GitInfo checkObj checkXdir remove_junk includes pre_commit

# Removes all compiled files: executables, objects, dependency files etc.
# Also deletes junk
clean:
	rm -f -v $(ALLEXES)
	rm -rf -v $(BD)/*.o $(BD)/*.d $(BD)/*.gc* $(BD)/*.txt $(BD)/*/
	rm -f -v *deleteme*

# Deletes junk outputs (e.g. '*deleteme') created by test suits
remove_junk:
	rm -f -v *deleteme*

# Make the 'ampsci.pdf' physics documentation (uses latex)
docs:
	( cd ./doc/tex && make )
	cp ./doc/tex/ampsci.pdf ./doc/ampsci.pdf
	( cd ./doc/tex && make clean)

# Make the doxygen code documentation (called 'docs' target)
doxy:
	doxygen ./doc/doxygen/Doxyfile
	( cd ./doc/latex && make )
	cp ./doc/latex/refman.pdf ./doc/documentation.pdf
	make docs
	cp ./doc/ampsci.pdf ./docs/ampsci.pdf 2>/dev/null || :
	( cd ./doc/latex && make clean)

# Extremely important
do_the_chicken_dance:
	@echo 'Why would I do that?'

# Runs clang format, and applies the changes to all c++ files
clang_format:
	@echo "Running clang format (whole project)"
	find $(SD)/ \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) ! -iname 'catch.hpp' | xargs clang-format -i -verbose


# Should be run before committing to dev
pre_commit:
	make includes
	make clang_format
	make ampsci tests
	./tests [unit]
	make remove_junk