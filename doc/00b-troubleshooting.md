\page troubleshooting Basic troubleshooting
\ingroup getting_started
\brief Basic troubleshooting and common errors

<!-- ## Common Compilation errors -->

A small collection of common errors users have seen, and how to beat them.
If you run into more exotic errors, let me know and I'll add them here.

### Makefile errors

If running 'make' doesn't compile the code, or you get a "nothing to be done" or an error like:

* ```make: *** No rule to make target '/main.o', needed by 'ampsci'.  Stop.```

you might be using an out-of-date Makefile.
Either re-run `configure.sh` or manually re-copy the updated makefile:

<div class="shell-block">
```shell
./configure.sh
```
</div>
or
<div class="shell-block">
```shell
cp ./doc/Makefile .
make
```
</div>

and try again.

### libgfortran error

Error linking to gfortran/lapack, may be something like this:

```text
ld: .../liblapack.a(xerbla.o): undefined reference to symbol '_gfortran_string_len_trim@@GFORTRAN_8'
ld: .../libgfortran.so.5: error adding symbols: DSO missing from command line
make: *** [build/buildTargets.mk:66: ampsci] Error 1
```

Sometimes occurs on certain systems (e.g., bunya).
In these cases, the fortran libraries need to be linked to explicitely.
This can be done by adding `-lgfortran` to the `LDLIBS` option in the Makefile

```Make
LDLIBS += -lgfortran
```

### OpenMP errors

* **error: unsupported option -fopenmp**
* **error: could not find <omp.h>**

* openmp (used for parallelisation) is not working. See [Compilation Details](\ref compilation) for some possible solutions.
* Quick fix: change '_OMPLIB=-fopenmp_' to '_OMPLIB=_' in Makefile, which turns off OpenMP

### GSL Errors

* **fatal error: gsl/<...>.h: No such file or directory** (or similar)
* **gsl** related linking/compilation error:

* Could not find required GSL libraries. Either they are not installed, or you need to link to them
* 1) Ensure GSL is installed (see above for instructions)
* 1) If GSL library is not installed in _/usr/local/_, you have to tell the compiler where to find the GSL files. Do this by setting the _GSL\_PATH_ option in Makefile. Common examples:
  * _GSL\_PATH=/opt/gsl/2.1/gnu_ # For UQ's getafix cluster
  * _GSL\_PATH=/usr/local/opt/gnu-scientific-library_ # For my macbook
  * Note: the exact path may differ for you, depending on where GSL was installed

* **error: too few arguments to function ‘int gsl_bspline_deriv_eval**

* This is because the code is linking to a _very_ old version of GSL. You might need to update GSL. If you have updated GSL (to at least version 2.0) and still get the message, the code is probably linking against the wrong version of GSL; see above to point the compiler to the correct version
* On HPC cluster, this might be because you haven't loaded the gsl module `module load gsl` (_after_ loading the gcc compiler)

### ld: Assertion failed: (resultIndex < sectData.atoms.size())

Full error message may look something like:

```text
0  0x102497648  __assert_rtn + 72
1  0x1023cbfac  ld::AtomPlacement::findAtom(unsigned char, unsigned long long, ld::AtomPlacement::AtomLoc const*&, long long&) const + 1204
2  0x1023e1924  ld::InputFiles::SliceParser::parseObjectFile(mach_o::Header const*) const + 15164
3  0x1023eee30  ld::InputFiles::parseAllFiles(void (ld::AtomFile const*) block_pointer)::$_7::operator()(unsigned long, ld::FileInfo const&) const + 420
4  0x185a37950  _dispatch_client_callout2 + 20
5  0x185a4aba0  _dispatch_apply_invoke + 176
6  0x185a37910  _dispatch_client_callout + 20
7  0x185a493cc  _dispatch_root_queue_drain + 864
8  0x185a49a04  _dispatch_worker_thread2 + 156
9  0x185be10d8  _pthread_wqthread + 228
ld: Assertion failed: (resultIndex < sectData.atoms.size()), function findAtom, file Relocations.cpp, line 1336.
collect2: error: ld returned 1 exit status
```

* This seems to be a bug with M1/M2 mac linker implementation with CommandLineTools version 15.0 (see [developer.apple.com/forums/thread/737707](https://developer.apple.com/forums/thread/737707))
* See also: [github.com/Homebrew/homebrew-core/issues/145991](https://github.com/Homebrew/homebrew-core/issues/145991)

* This issue seems to have been resolved in version 15.1 - so run an update if you can.

* If you can't update, it seems this can be fixed by adding the following to the Makefile:

```Make
LARGS=-Wl,-ld_classic
```

You shouldn't need to make clean first, but if it doesn't work, try that too.

### Others

* Sometimes, the compiler will not be able to find the correct libraries (particular, e.g., on clusters). In this case, there are two options in the Makfefile: **CARGS ?=** and **LARGS ?=**
