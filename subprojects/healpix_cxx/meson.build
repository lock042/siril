################################################################################
# Meson main file
################################################################################

################################################################################
# Define project

project('healpix_cxx', 
  ['cpp'],
  version : '8.3.0',
  meson_version: '>=0.53.0',
  license : 'GPL-2.0-only',
  default_options : [
    'c_std=gnu11',
    'cpp_std=c++11'
    ],
)

################################################################################
# Meson modules and compilers

## Compilers
cpp = meson.get_compiler('cpp')

################################################################################
# Configuration

## General and package configuration
conf_data = configuration_data()

################################################################################
# Generate configuration files

################################################################################
# List sources

inc = include_directories('.')

src_files = ['error_handling.cc',
            'geom_utils.cc',
            'healpix_base.cc',
            'healpix_tables.cc',
            'pointing.cc',
            'string_utils.cc',
]

################################################################################
# Build healpix_cxx_lib library

healpix_cxx_lib = static_library('healpix_cxx_lib', src_files, include_directories: '.')

################################################################################
# Dependency wcslib

healpix_cxx_dep = declare_dependency(include_directories: '.', link_with : healpix_cxx_lib)
