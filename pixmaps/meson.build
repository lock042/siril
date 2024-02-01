################################################################################
# Meson file for pictures/icons/logo
################################################################################

################################################################################
# List icons

siril_icons = [
  'siril.svg'
]

################################################################################
# Install

install_subdir('icons/hicolor', install_dir : siril_share + '/icons')
install_data(siril_icons,
             install_dir : siril_datadir + '/pixmaps/' )
install_subdir('plot_background', install_dir : siril_datadir + '/pixmaps/')
