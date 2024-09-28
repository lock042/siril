This directory contains tests for siril commands written in python to ease the
process, using pysiril that connects to siril-cli.

The criterion C tests can be used for unit testing, those will be more for
command testing and non regression checks in important parts of large pieces of
code.

Some lightweight images can be put there too.

Files and test function names must start with 'test_'.

Install pysiril (see https://siril.org/tutorials/pysiril/), install pytest. To
run the tests, run the command pytest in this directory.

TODO:
* configurable siril-cli path (now works only with autotools)
* only one siril-cli instance that doesn't use the main config file to not mess
  with it and use on defaults on every start
* write many tests
* add to the CI?
