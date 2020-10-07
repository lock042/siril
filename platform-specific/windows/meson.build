################################################################################
# Meson file
################################################################################

################################################################################
# Compile resource file

gimp_app_rc = configure_file(
  input : 'siril.rc.in',
  output: 'siril.rc',
  configuration: versionconfig,
)

win_res += winm.compile_resources('siril.rc', depend_files: ['siril.ico', 'text-x-seq.ico'])