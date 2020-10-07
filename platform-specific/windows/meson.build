################################################################################
# Meson file
################################################################################

################################################################################
# Compile resource file

siril_rc = configure_file(
  input : 'siril.rc.in',
  output: 'siril.rc',
  configuration: versionconfig,
)

inno_setup = configure_file(
  input : 'WindowsInnoSetup.iss.in',
  output : 'WindowsInnoSetup.iss',
  configuration : conf_data,
)

win_res += winm.compile_resources(
   siril_rc,
   args: [
      '--define', 'ORIGINALFILENAME_STR="@0@"'.format(meson.project_name()+'.exe'),
      '--define', 'INTERNALNAME_STR="@0@"'    .format(meson.project_name()),
      '--define', 'TOP_SRCDIR="@0@"'          .format(meson.source_root()),
    ],
   depend_files: ['/platform-specific/windows/siril.ico', '/platform-specific/windows/text-x-seq.ico'],
)