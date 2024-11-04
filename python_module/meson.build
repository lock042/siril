# List of Python files to install
python_files = [
  'setup.py',
  'siril/__init__.py',
  'siril/connection.py',
  'siril/exceptions.py',
  'siril/models.py',
  'siril/shm.py',
  'siril/translations.py',
  # Add other .py files here
]

# List of markdown files
markdown_files = [
  'README.md',
  'locale/README.md',
  # Add other .md files here
]

# List of message catalog files
#message_files = [
#  'locale/fr/LC_MESSAGES/siril.mo',
#  'locale/es/LC_MESSAGES/siril.mo',
  # Add other .mo files here
#]

# Install all files
install_data(python_files + markdown_files, # + message_files,
            install_dir: siril_datadir + '/python_module')

# Ensure directories are created
install_subdir('siril',
  install_dir: siril_datadir + '/python_module',
  exclude_files: ['*.py', '*.md', '*.mo']
)

install_subdir('locale',
  install_dir: siril_datadir + '/python_module',
  exclude_files: ['*.py', '*.md', '*.mo']
)
