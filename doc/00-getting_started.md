\page getting_started Getting started

# Getting the code

* The code is available on GitHub: [github.com/benroberts999/ampsci](https://github.com/benroberts999/ampsci).
* It's recommended to use `git` to grab the code, since then you can most easily pull updates (which happen regularly). This is already installed on most linux distributions, and comes with homebrew on mac, but may need to be installed, e.g.:

<div class="shell-block">
```shell
  sudo apt install git
```
</div>

* If you have GitHub ssh keys setup (see [github docs for ssh keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)), get the code using:

<div class="shell-block">
```shell
  git clone git@github.com:benroberts999/ampsci.git
```
</div>

* Or, if you don't have them set up:

<div class="shell-block">
```shell
  git clone https://github.com/benroberts999/ampsci.git
```
</div>

* With this, you can regularly pull any updates using pull (remember to re-comile!):

<div class="shell-block">
```bash
git pull
```
</div>

* Alternatively, you may just download the source code [ampsci.zip](https://github.com/benroberts999/ampsci/archive/refs/heads/main.zip) [not recommended]

## Quick-start: compiling ampsci

The following simple method should work for most cases. More detailed compilation instructions/options are provided: [Compilation Details](\ref compilation)

* Easiest method to quickly get running is to use provided shell script `install-dependencies.sh`, which installs the required packages/compilers
  * This might not work on all systems, meaning a manual setup will be required.
  * If on linux, you will need sudo privileges to run `install-dependencies.sh`
  * On Mac - the easiest method is to use homebrew [https://brew.sh/](https://brew.sh/) to install the required packages (`install-dependencies.sh` will do this for you)

<div class="shell-block">
```shell
./install-dependencies.sh
```
</div>

* Alternatively, directly install the packages. E.g., on Ubuntu:

<div class="shell-block">
```shell
  sudo apt install g++ liblapack-dev libblas-dev libgsl-dev make
```
</div>

* Or on Mac (once [homebrew](https://brew.sh/) is installed):

<div class="shell-block">
```shell
  brew install gcc gsl
```
</div>

* Then, run `setup.sh` which will set up the Makefile and do the initial compilation

<div class="shell-block">
```shell
./setup.sh
```
</div>

* These scripts will setup the Makefile and compile ampsci. After these have run, they do not need to be run again. You can re-compile ampsci after this using the makefil simply by running `make` from the command line. You will need to re-run make if you make make any changes to the code, or pull changes down from GitHub:

<div class="shell-block">
```shell
make
```
</div>

* The provided file _Makefile_ has some basic compilation options. It's currently set up to work on most linux systems; you may need to change a few options for others (see detailed [Compilation Details](\ref compilation))

* After that, you can run the first example input file:

<div class="shell-block">
```shell
  cp ./doc/examples/ampsci.in ./
  ./ampsci ampsci.in
```
</div>

## Documentation

* The documentation is provided on this website
* A pdf with detailed physics descriptions of the methods is also provided: [ampsci.pdf](https://ampsci.dev/ampsci.pdf)
* If you want, you can build the html and pdf documentaion locally using the Makefile:
  * This requires doxygen for the html docs (`apt install doxygen`), and pdflatex for the pdf docs.

<div class="shell-block">
```shell
make docs
```
</div>

-----

## See also

* \subpage compilation - \copybrief compilation

* \subpage troubleshooting - \copybrief troubleshooting

* \subpage examples - \copybrief examples
