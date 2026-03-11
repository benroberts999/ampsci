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

.PHONY: all full clean remove_junk docs clang_format pre_commit includes license-year force do_the_chicken_dance

clean:
	rm -fv $(ALLEXES)
	rm -rfv $(BUILD)/* $(BUILD)/*/

## Remove unit test outputs
remove_junk:
	rm -f -v *deleteme*

## Builds the ampsci.pdf docs using tex, and the html using doxygen, if available
docs:
	( cd ./doc/tex && $(MAKE) 2>/dev/null || : )
	cp ./doc/tex/ampsci.pdf ./doc/ampsci.pdf 2>/dev/null || :
	cp ./LICENSE ./doc/LICENSE.md
	doxygen ./doc/doxygen/Doxyfile 2>/dev/null && \
		cp ./doc/doxygen/ampsci.html ./doc/ampsci-documentation.html 2>/dev/null || :
	cp ./doc/tex/ampsci.pdf ./doc/html/ampsci.pdf 2>/dev/null || :
	cp -r ./doc/img ./doc/html/ 2>/dev/null || true
	cp ./doc/doxygen/search.js ./doc/html/search/search.js
	echo "<script id=\"searchdata\" type=\"text/xmldata\">" >> ./doc/html/search.html
	cat ./doc/searchdata.xml >> ./doc/html/search.html
	echo "</script>" >> ./doc/html/search.html

do_the_chicken_dance:
	@echo 'Why would I do that?'

## Run clang format on entire project
## External libraries (fmt, catch2, json) are excluded via their own .clang-format files
## Force specific clang-format version to avoid many small diffs
CLANG_FORMAT := clang-format-14
clang_format:
	@echo "Running clang-format (whole project)"
	@find $(SRC) \
	  -type f \( -iname '*.cpp' -o -iname '*.hpp' -o -iname '*.ipp' \) -print \
	  | xargs -r $(CLANG_FORMAT) -i -verbose

## Makes the include.hpp files, and runs clang format on them
includes:
	$(MAKE) license-year
	./$(SRC)/tools/build_includes.sh
	@echo "Running clang format (on includes)"
	find $(SRC)/ -iname 'include.hpp' | xargs clang-format -i -verbose

## Updates Licence year
license-year:
	@year=$$(date +%Y); \
	sed -i.bak -E "s/(Copyright \(c\) )([0-9]{4})(-[0-9]{4})?/\1\2-$${year}/" LICENSE; \
	rm -f LICENSE.bak; \
	echo "Updated license year to $${year}"

## All things should be done before major commits
pre_commit:
	$(MAKE) includes
	$(MAKE) clang_format
	$(MAKE) license-year
	$(MAKE) ampsci tests
	./test [unit]
	rm -f -v *deleteme*
