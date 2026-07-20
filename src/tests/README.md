# SIRIL test suites

## Prerequisites

* Siril test suites use [criterion](criterion https://criterion.readthedocs.io/en/master/intro.html). If your distribution does not have any packages, you need to compile it by following instruction [here](https://github.com/Snaipe/Criterion.git).

* Build Siril:

    mkdir _build
    meson --buildtype debug _build -Dcriterion=true
    ninja -C _build

## Running tests

### What is inside the test directory?
There are different kinds of files in this directory:
- some unit tests, named `something_test.c`,
- a performance test for sorting algorithms for median computation in
  `sorting.c`,
- a program that can be used to compare FITS files, to verify that an algorithm
  always computes the same thing for example, in the file `compare_fits.c` 
- a build script `script.sh` for the two that are not unit tests.

The build script and the two executables existed before the meson build system,
but the meson build has not been updated to generate them. Since they depend on
siril's code and we don't want to pull all the files here, we had to redefine
the functions used in the direct dependencies in dummy.c, they should not be
used in the tests.

### Using criterion for unit tests

Criterion is integrated to meson and it's the right way to run the unit tests:

    # run all tests
    meson test

    # run the test ser_test and show the output
    meson test -v ser_test

### Using test script

To compile the test programs, compile siril with the autotools/make method,
then run ./build.sh in src/tests.
Since sorting makes some performance tests, siril has to be compiled with -O2
to have realistic values.

If build error occurs, check that the basic build script has all required
package links for your OS and options for your compiler.
If link error occurs, add the missing functions in dummy.c

This is very rarely used, it will probably not work.

## Debugging tests

The script creates executables for some tests, which can be debugged like any other.
For the meson and criterion method, there are a few useful commands:

    # run the test 'ser_test' with valgrind and print logs to stdout
    meson test --wrap=valgrind ser_test --print-errorlogs

    # run the test ser_test with gdb (just run it on prompt)
    meson test --gdb ser_test

Unfortunately meson --gdb doesn't seem to work, it never breaks.
The criterion debug method with gdbserver doesn't work either.
For this reason, some tests can also be compiled as a regular executable, with
the criterion calls replaced by a macro that displays errors. But as for the
previous section, the build script is rarely tested.

