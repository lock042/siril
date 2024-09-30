# SIRIL test suites

## Prerequisites

* Siril test suites use [criterion](criterion
  https://criterion.readthedocs.io/en/master/intro.html), install it. If your
distribution does not have any packages, you need to compile it by following
instruction [here](https://github.com/Snaipe/Criterion.git).

* The automated test system is only available through meson, Build with
  criterion enabled by passing `-Dcriterion=true`, it's disabled by default:
    mkdir _build
    meson _build -Dcriterion=true
    ninja -C _build
    ninja -C _build install

* Run the tests with:
    meson test -v

## Running tests

As shown above, the integrated unit test system is available through meson. For
different settings, see this [meson documentation
page](https://mesonbuild.com/Unit-tests.html). Tests can be run several times,
selected to run only by their name, run wrapped in valgrind, and so on.

Originally, the tests were compiled as executables using autotools and the
build.sh script in this directory. There are also two special files in here that
are not unit tests:
- `compare_fits` is a program that can be used to compare FITS files, to verify
  that an algorithm always computes the same thing for example
- `sorting` is a unit test on the three sorting implementations that provide the
  median, but it also contains a performance evaluation between them, so it
takes more time than a simple unit test.

### Using test script (autotools)

Because on some OS the debugging of tests in meson is not supported, tests
can also be compiled as a regular executable, with the criterion calls replaced
by a macro that displays errors.

Since the executables depend on siril's code and we don't want to pull all the
files here, we had to redefine the functions used in the direct dependencies in
dummy.c, they should not be used in the tests. This is rarely updated and
probably doesn't work as is.

To compile the test programs, compile siril with the autotools/make method,
then run ./build.sh in src/tests.
Since sorting makes some performance tests, siril has to be compiled with -O2
to have realistic values.

If build error occurs, check that the basic build script has all required
package links for your OS and options for your compiler.
If link error occurs, or a segmentation fault happens on launch, add the missing
functions in dummy.c

## Debugging scripts

With the autotools way, the script creates executables for the tests which can
be debugged like any program. With the meson and criterion method, there are a
few useful commands:

    # run the test 'ser_test' with valgrind and print logs to stdout
    meson test -C _build --wrap=valgrind ser_test --print-errorlogs

Unfortunately meson --gdb doesn't work, it never breaks, even if there is a
segmentation fault in the test. Debugging with gdb using the method [documented
in criterion](https://criterion.readthedocs.io/en/master/debug.html), where the
test is started in gdb but another remote gdb session is required to access it,
works on linux and probably mac.

