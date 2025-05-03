# List of Python files to install (outside of siril subdirectory)
python_files = [
  'pyproject.toml',
  'MANIFEST.in',
]

# List of markdown files to install
markdown_files = [
  'README.md',
  'locale/README.md',
]

# Install individual files (not inside subdirectories)
install_data(python_files + markdown_files,
            install_dir: siril_datadir + '/python_module')

# Install the entire 'sirilpy' subdirectory (this will handle all .py files in sirilpy/)
install_subdir('sirilpy',
  install_dir: siril_datadir + '/python_module'
)

# Install the 'locale' directory, excluding Python and markdown files
install_subdir('locale',
  install_dir: siril_datadir + '/python_module',
  exclude_files: ['*.py', '*.md']
)
