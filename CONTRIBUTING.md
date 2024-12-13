# Contributing to Siril

First, thanks for taking the time to contribute!

## How Can I Contribute?

Siril is Free Software and you are welcome to contribute to this project. There are many ways to do it:

* develop new features,
* report bugs (errors in the program),
* test existing features and provide feedback,
* add or improve the documentation or tutorials,
* translate Siril to your own language,
* translate the documentation,
* donate

### Getting last version

In order to try the last development version we recommend to compile the sources (see [README](README.md)). Nevertheless it is possible to get current builds for all platforms directly from our CI by following these direct links:

* GNU/Linux ([x86_64](https://gitlab.com/free-astro/siril/-/jobs/artifacts/master/download?job=appimage))
* Windows cross-build ([x86_64](https://gitlab.com/free-astro/siril/-/jobs/artifacts/master/download?job=win64))
* macOS ([arm64](https://gitlab.com/free-astro/siril/-/jobs/artifacts/master/download?job=siril-macos:%20[macosarm]), [x86_64](https://gitlab.com/free-astro/siril/-/jobs/artifacts/master/download?job=siril-macos:%20[shared-macos-amd64]))

 **Test builds are for testing purpose only. They have not been human-tested, it relies on regularly modified development code. So please do not use it for production!**

### Reporting Bugs

Reporting the bugs that you will encounter is very important to the development, it helps the developers to make Siril more stable and more bug free. If you have some programming skills you can attach a patch to your bug report, we will be happy to apply it.

#### Before Submitting A Bug Report

Before creating bug reports, please check [this list](https://gitlab.com/free-astro/siril/issues) as you might find out that you don't need to create one. Also, make sure you are using last stable version or last git version before submitting the ticket.
When you are creating a bug report, please include as many details as possible. Fill out [the required template](https://gitlab.com/free-astro/siril/blob/master/.gitlab/issue_templates/bug.md), the information it asks for helps us resolve issues faster.

### Getting Your Code Changes Merged

If you would like to fix any bug or add a new feature, you can refer to our [contributor guidelines](https://siril-contrib-doc.readthedocs.io/en/latest/).
In particular, you will need to have a look at the section about [submitting your work](https://siril-contrib-doc.readthedocs.io/en/latest/SubmitingWork.html).

### Translation

If you're interested in contributing to the translation of the application and documentation, 
we encourage you to use [Weblate](https://weblate.siril.org/projects/siril/). Weblate is a 
powerful web-based translation tool that allows for easy collaboration and efficient 
translation workflow.

### Formatting

* Formatting code

Code formatting follows [K&R style](https://en.wikipedia.org/wiki/Indentation_style#K&R_style).

* Formatting commands help

When adding commands, you will need to add a description of what it does and which arguments it takes, if any, in `src/core/command_def.h`. All arguments names must be put in bold and possible string values between quotes. The last sentence of the command help must not be ended by a period.

### Donate

If you like the software, please help us by contributing with the [Donate button](https://www.siril.org/#support-us) of the website. Siril takes us a lot of time and we still have to pay for the servers.
