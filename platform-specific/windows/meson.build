################################################################################
# Meson file
################################################################################

################################################################################

siril_rc = configure_file(
  input : 'siril.rc.in',
  output: 'siril.rc',
  configuration: conf_data,
)


# Compile resource file

win_res += winm.compile_resources(
    siril_rc,
    args: [
      '--define', 'ORIGINALFILENAME_STR="@0@"'.format(meson.project_name()+'.exe'),
      '--define', 'INTERNALNAME_STR="@0@"'    .format(meson.project_name()),
      '--define', 'TOP_SRCDIR="@0@"'          .format(meson.project_source_root()),
    ],
    include_directories: [
      rootInclude, srcInclude,
    ],
		depend_files: ['siril.ico', 'text-x-seq.ico']
	)