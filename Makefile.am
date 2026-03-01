## Process this file with automake to produce Makefile.in

export CC CFLAGS CPPFLAGS

SUBDIRS = subprojects \
	po \
	src \
	data \
	pixmaps \
	platform-specific \
	python_module

# TODO add a conditionnal and a configure switch
# to be able to choose the language. Also, this should
# be a recussive make with a separate file list.
scriptsdir = $(pkgdatadir)/scripts
scripts_DATA = scripts/OSC_Preprocessing.ssf \
    scripts/Mono_Preprocessing.ssf \
	scripts/OSC_Extract_Ha.ssf \
	scripts/OSC_Extract_HaOIII.ssf \
	scripts/OSC_Preprocessing_BayerDrizzle.ssf

doc_DATA = AUTHORS \
	LICENSE.md \
	3rdparty_licenses/LICENSE_sleef.txt \
	3rdparty_licenses/LICENSE_zlib.txt \
	3rdparty_licenses/GPL-2.0-or-later.txt

man_MANS = siril.1

if !PLATFORM_WIN32
install-exec-hook:
	(cd $(DESTDIR)$(mandir)/man1; \
	rm -f siril-cli.1; \
	$(LN_S) siril.1 siril-cli.1)

uninstall-hook:
	(cd $(DESTDIR)$(mandir)/man1; rm -f siril-cli.1)
endif

ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = config.rpath \
             git-version.h \
             siril_resource.h
AUTOMAKE_OPTIONS = foreign

generated_sources = \
	git-version.h \
	siril_resource.h

# Build git-version.h before anything in the subdirs as this is needed
# in the about and debug dialog (app/) and in the debug tool (tools/).
BUILT_SOURCES = $(generated_sources)

git-version.h: update-git-version-header
	@if test -e "$(top_srcdir)/.git"; then \
	  git_version="`git --git-dir=$(top_srcdir)/.git describe --always`"; \
	  git_version_abbrev="`git --git-dir=$(top_srcdir)/.git rev-parse --short HEAD`"; \
	  git_last_commit_year="`git --git-dir=$(top_srcdir)/.git log -n1 --reverse --pretty=%ci | cut -b 1-4`"; \
	elif test ! -f "$@"; then \
	  git_version="Unknown"; \
	  git_version_abbrev="$$git_version"; \
	  git_last_commit_timestamp=-1; \
	  git_last_commit_year="`date -u '+%Y'`"; \
	fi; \
	if test -n "$$git_version"; then \
	  echo "#ifndef __GIT_VERSION_H__" > "$@.tmp"; \
	  echo "#define __GIT_VERSION_H__" >> "$@.tmp"; \
	  echo "#define SIRIL_GIT_VERSION_ABBREV \"$$git_version_abbrev\"" >> "$@.tmp"; \
	  echo "#define SIRIL_GIT_LAST_COMMIT_YEAR \"$$git_last_commit_year\"" >> "$@.tmp"; \
	  echo "#endif /* __GIT_VERSION_H__ */" >> "$@.tmp"; \
	fi
	@if ( test -f "$@.tmp" && test -f "$@" && cmp "$@.tmp" "$@" > /dev/null ); then \
	  rm -f "$@.tmp"; \
	elif test -f "$@.tmp"; then \
	  mv "$@.tmp" "$@"; \
	  echo "  git HEAD changed: $@ regenerated"; \
	fi

siril_resource.h: src/siril_resource.xml
	$(AM_V_GEN)
	glib-compile-resources \
	--sourcedir=$(top_builddir) \
	--sourcedir=src \
	--target=src/$@ \
	--generate-header \
	--c-name siril_resource \
	$(srcdir)/src/siril_resource.xml

.PHONY: update-git-version-header siril_resource.h
