# List of Python files to install (outside of siril subdirectory)
python_files = [
  'pyproject.toml',
  'MANIFEST.in',
]

# Install individual files (not inside subdirectories)
install_data(python_files,
            install_dir: siril_datadir + '/python_module'
)
            
install_data(['README.md'],
            install_dir: siril_share + 'doc/siril/python_module'
)

# Create the locale directory in doc and install its README.md
install_data(['locale/README.md'],
            install_dir: siril_share + 'doc/siril/python_module/locale'
)

# Install the entire 'sirilpy' subdirectory (this will handle all .py files in sirilpy/)
install_subdir('sirilpy',
            install_dir: siril_datadir + '/python_module'
)