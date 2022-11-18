# All makefile targets

################################################################################
#Allow exectuables to be placed in another directory:
DEFAULTEXES = $(addprefix $(XD)/, \
 ampsci \
)

ALLEXES = $(addprefix $(XD)/, \
 ampsci dmeXSection tests \
)

CHECKS = GitInfo checkObj checkXdir

#Default make rule:
all: $(CHECKS) $(DEFAULTEXES)
full: $(CHECKS) $(ALLEXES)

################################################################################
# Automatically generate dependency files for each cpp file, + compile:
# Put built objects in subdirectories based on build mode and compiler

SUBBD=$(Build)/$(word 1,$(CXX))

# Auto rule for all cpp files
$(BD)/$(SUBBD)/%.o: $(SD)/%.cpp
	@mkdir -p $(@D)
	$(COMP)

# Force version.o to build each time: ensure updated git+version info!
# Not the best solution, but works?
.PHONY: force
$(BD)/$(SUBBD)/version/version.o: $(SD)/version/version.cpp force
	@mkdir -p $(@D)
	$(COMP) $(GITFLAGS)

# include the dependency files
-include $(BD)/$(SUBBD)/*.d $(BD)/$(SUBBD)/*/*.d

################################################################################
# List all objects in sub-directories (i.e., that don't conatin a main())
# e.g., for each src/sub_dir/file.cpp -> build/sub_dir/file.o

# Each test file (except main()):
TEST_SRC_FILES := $(wildcard $(SD)/*/*.tests.cpp)
TEST_OBJS := $(subst $(SD),$(BD)/$(SUBBD),$(subst .cpp,.o,$(TEST_SRC_FILES)))

#Each non-test source file (except main())
SRC_FILES := $(filter-out $(TEST_SRC_FILES), $(wildcard $(SD)/*/*.cpp))
OBJS := $(subst $(SD),$(BD)/$(SUBBD),$(subst .cpp,.o,$(SRC_FILES)))

################################################################################
# Link + build all final programs

$(XD)/ampsci: $(BD)/$(SUBBD)/ampsci.o $(OBJS)
	$(LINK)

$(XD)/tests: $(OBJS) $(TEST_OBJS)
	$(LINK)

$(XD)/dmeXSection: $(BD)/$(SUBBD)/dmeXSection.o $(OBJS)
	$(LINK)


################################################################################
# Add git version info to compile flags
NOW:=$(shell date +%Y-%m-%d' '%H:%M' '%Z 2>/dev/null)
GITFLAGS=-D GITREVISION="$(shell git rev-parse --short HEAD 2>/dev/null)"
GITFLAGS+=-D GITBRANCH="$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)"
GITFLAGS+=-D GITMODIFIED="$(shell git status -s 2>/dev/null)"
GITFLAGS+=-D CXXVERSION="$(shell $(CXX) --version 2>/dev/null | sed -n '1p' | sed s/'('/'['/ | sed s/')'/']'/)"
GITFLAGS+=-D COMPTIME="$(NOW)"

################################################################################
################################################################################

GitInfo:
	@echo Git Branch: $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)
	@echo Git Revision: $(shell git rev-parse --short HEAD 2>/dev/null)
	@echo Modified: $(shell git status -s 2>/dev/null)

checkObj:
	@if [ ! -d $(BD) ]; then \
	  echo '\n ERROR: Directory: '$(BD)' doesnt exist - please create it!\n'; \
	  false; \
	fi

checkXdir:
	@if [ ! -d $(XD) ]; then \
		echo '\n ERROR: Directory: '$(XD)' doesnt exist - please create it!\n'; \
		false; \
	fi

.PHONY: clean docs doxy do_the_chicken_dance GitInfo checkObj checkXdir remove_deleteme
clean:
	rm -f -v $(ALLEXES)
	rm -rf -v $(BD)/*.o $(BD)/*.d $(BD)/*.gc* $(BD)/*/
	rm -f -v *deleteme*
# Make the 'ampsci.pdf' physics documentation
docs:
	( cd ./doc/tex && make )
	cp ./doc/tex/ampsci.pdf ./doc/ampsci.pdf
	( cd ./doc/tex && make clean)
# Make the doxygen code documentation
doxy:
	doxygen ./doc/doxygen/Doxyfile
	( cd ./doc/latex && make )
	cp ./doc/latex/refman.pdf ./doc/documentation.pdf
	make docs
	cp ./doc/ampsci.pdf ./docs/ampsci.pdf 2>/dev/null || :
	( cd ./doc/latex && make clean)
do_the_chicken_dance:
	@echo 'Why would I do that?'
clang_format:
	clang-format -i src/*.cpp src/*/*.cpp src/*/*.hpp
remove_deleteme:
	rm -f -v *deleteme*