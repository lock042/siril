################################################################################
# Meson file
################################################################################

################################################################################
# Define platform specific variables

win_res = []

################################################################################
# Then add the correct folder

if platform_linux or platform_bsd
  subdir('linux')
elif platform_windows
  subdir('windows')
elif platform_osx
  subdir('macos')
else
  error('Unsupported OS.')
endif