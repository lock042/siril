#SIRIL test suites
=====

##Prerequisites
-------
* First, you need to compile Siril

    mkdir _build
    meson --buildtype release _build
    cd _build
    ninja
    ninja install

* Then, Siril test suites use [criterion](criterion https://criterion.readthedocs.io/en/master/intro.html). If your distribution does not have any packages, you need to compile it by following instruction [here](https://github.com/Snaipe/Criterion.git).
    


##Running tests
------
###What are inside test directory?
There are different kinds of files in this directory:
- `compare_fits` is a program that can be used to compare FITS files, to verify
  that an algorithm always computes the same thing for example
- `sorting` is a unit test on the three sorting implementations that provide the
  median. It also contains a performance evaluation between them.

Other files are used for the build of these executables. Since they depend on
siril's code and we don't want to pull all the files here, we had to redefine
the functions used in the direct dependencies in dummy.c, they should not be
used in the tests.

###Using criterion

    meson test (or meson test -C _build from src)
    meson test -v

###Using test script

To compile the test programs, compile siril with the autotools/make method,
then run ./build.sh in src/tests.
Since sorting makes some performance tests, siril has to be compiled with -O2
to have realistic values.

If build error occurs, check that the basic build script has all required
package links for your OS and options for your compiler.
If link error occurs, add the missing functions in dummy.c

## Debugging scripts

The script creates executables for some tests, which can be debugged like any other.
For the meson and criterion method, there are a few useful commands:

    # run the test 'ser_test' with valgrind and print logs to stdout
    meson test --wrap=valgrind ser_test --print-errorlogs

    # run the test ser_test with gdb (just run it on prompt)
    meson test --gdb ser_test

Unfortunately meson --gdb doesn't work on raspbian, it never breaks.
The criterion debug method with gdbserver doesn't work either.
For this reason, some tests can also be compiled as a regular executable, with
the criterion calls replaced by a macro that displays errors.

