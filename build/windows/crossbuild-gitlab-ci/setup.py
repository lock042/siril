#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This file is part of crossroad.
# Copyright (C) 2013-2014 Jehan <jehan at girinstud.io>
#
# crossroad is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# crossroad is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with crossroad.  If not, see <http://www.gnu.org/licenses/>.

import distutils.command.build
import distutils.command.build_scripts
import distutils.command.install_data
import distutils.command.install_scripts
import gzip
import sys
import os
import stat
import subprocess
import shutil
import configparser

version = '0.9.0'
deactivated_platforms = ['arm', 'arm-gnu']
srcdir   = os.path.dirname(os.path.realpath(sys.argv[0]))
builddir = os.getcwd()

use_setuptools = False

# if 'USE_SETUPTOOLS' in os.environ or 'setuptools' in sys.modules:
#     # I don't like to unreference modules out of their namespaces.
#     # Unfortunately it seems pip would need setuptools instead of the
#     # core distutils. Thus I import setup and install from setuptools
#     # when it is requested or already loaded (ex: in pip).
#     use_setuptools = True
#     try:
#         from setuptools import setup
#         from setuptools.command.install import install
#     except ImportError:
#         use_setuptools = False
print('Expe setup.py')

if not use_setuptools:
    from distutils.core import setup
    from distutils.command.install import install

class build_man(distutils.core.Command):
    '''
    Build the man page.
    '''

    description = 'build the man page'
    user_options = []

    def run(self):
        self.check_dep()
        self.create_build_tree()
        self.build()
        self.compress_man()
        self.clean()

    def initialize_options(self):
        self.cwd = None

    def finalize_options(self):
        self.cwd = os.getcwd()

    def check_dep(self):
        '''
        Check build dependencies.
        '''
        if shutil.which('rst2man') is None:
            sys.stderr.write('`rst2man` is a mandatory building dependency. '
                             'You will probably find it either in a '
                             '`python2-docutils` or `python3-docutils` package.\n')
            sys.exit(os.EX_CANTCREAT)

    def create_build_tree(self):
        '''
        Create a build tree.
        '''
        try:
            os.makedirs('build/share/man/man1', exist_ok=True)
            os.makedirs('build/doc', exist_ok=True)
        except os.error:
            sys.stderr.write('Build error: failure to create the build/ tree. Please check your permissions.\n')
            sys.exit(os.EX_CANTCREAT)

    def build(self):
        '''
        Create the manual.
        '''
        try:
            shutil.copyfile(os.path.join(srcdir, 'doc/crossroad.rst'),
                            os.path.join(builddir, 'build/doc/crossroad.rst'))
            update_scripts('build/doc')
            subprocess.check_call(["rst2man", "build/doc/crossroad.rst",
                                   "build/share/man/man1/crossroad.1"])
        except subprocess.CalledProcessError:
            sys.stderr.write('Build error: `rst2man` failed to build the man page.\n')
            sys.exit(os.EX_CANTCREAT)

    def compress_man(self):
        '''
        Compress the man.
        '''
        with open('build/share/man/man1/crossroad.1', 'rb') as manual:
            with gzip.open('build/share/man/man1/crossroad.1.gz', 'wb') as compressed:
                compressed.writelines(manual)

    def clean(self):
        os.unlink('build/share/man/man1/crossroad.1')

class my_build(distutils.command.build.build):
    '''
    Override the build to have some additional pre-processing.
    '''

    def run(self):
        # Add manual generation in build.
        self.run_command('man')
        # Move files without modification at build time
        # This allows renaming mostly.
        try:
            os.makedirs('build/bin', exist_ok = True)
            os.makedirs('build/platforms', exist_ok = True)
            os.makedirs('build/share/crossroad/scripts/shells/bash', exist_ok = True)
            os.makedirs('build/share/crossroad/scripts/cmake', exist_ok = True)
            os.makedirs('build/share/crossroad/scripts/meson', exist_ok = True)
        except os.error:
            sys.stderr.write('Build error: failure to create the build/ tree. Please check your permissions.\n')
            sys.exit(os.EX_CANTCREAT)
        shutil.copyfile(os.path.join(srcdir, 'src/crossroad.py'), 'build/bin/crossroad')
        shutil.copyfile(os.path.join(srcdir, 'src/in-crossroad.py'), 'build/share/crossroad/scripts/in-crossroad.py')
        shutil.copy(os.path.join(srcdir, 'scripts/shells/environment.sh'), 'build/share/crossroad/scripts/shells/')
        for f in os.listdir(os.path.join(srcdir, 'platforms/env/')):
            if f[-5:] == '.conf' and f[:-5] not in deactivated_platforms:
                config = configparser.ConfigParser()
                config.optionxform = str
                config.read([os.path.join(srcdir, 'platforms/env', f)])
                if not config.has_section('Platform') or not config.has_option('Platform', 'shortname') or \
                   not config.has_option('Platform', 'nicename') or not config.has_option('Platform', 'host'):
                    sys.stderr.write('Build error: file {} miss required options.\n'. f)
                    sys.exit(os.EX_CANTCREAT)
                shortname = config.get('Platform', 'shortname')
                nicename = config.get('Platform', 'nicename')
                sys.stdout.write('Configuring platform "{}"\n'.format(nicename))

                env_variables = '\nexport CROSSROAD_PLATFORM="{}"\n'.format(config.get('Platform', 'shortname'))
                env_variables += 'export CROSSROAD_PLATFORM_NICENAME="{}"\n'.format(config.get('Platform', 'nicename'))
                env_variables += 'export CROSSROAD_HOST="{}"\n'.format(config.get('Platform', 'host'))

                if config.has_option('Platform', 'word-size'):
                    env_variables += 'export CROSSROAD_WORD_SIZE="{}"\n'.format(config.getint('Platform', 'word-size'))
                if config.has_section('Environment'):
                    custom_env_vars = config.items('Environment')
                    for (env_var, env_val) in custom_env_vars:
                        env_variables += 'export {}="{}"\n'.format(env_var, env_val)

                # Platform file.
                shutil.copy(os.path.join(srcdir, 'platforms/modules/', shortname + '.py'), 'build/platforms')
                if shortname != 'native':
                    # Cmake file.
                    shutil.copy(os.path.join(srcdir, 'platforms/cmake/', 'toolchain-' + shortname + '.cmake'), 'build/share/crossroad/scripts/cmake/')
                    # Meson cross build definition files.
                    shutil.copy(os.path.join(srcdir, 'platforms/meson/', 'toolchain-' + shortname + '.meson'), 'build/share/crossroad/scripts/meson/')
                # Bash startup file.
                built_bashrc = os.path.join('build/share/crossroad/scripts/shells/bash/', 'bashrc.' + shortname)
                shutil.copyfile(os.path.join(srcdir, 'scripts/shells/bash/bashrc.template'), built_bashrc)
                try:
                    fd = open(built_bashrc, 'r+')
                    contents = fd.read()
                    contents = contents.replace('@ENV_VARIABLES@', env_variables)
                    fd.truncate(0)
                    fd.seek(0)
                    fd.write(contents)
                    fd.flush()
                    fd.close()
                except IOError:
                    sys.stderr.write('"{}" failed to update. Check your permissions.\n'.format(built_bashrc))
                    sys.exit(os.EX_CANTCREAT)
                # Zsh startup files.
                zsh_startup_dir = 'build/share/crossroad/scripts/shells/zsh.' + shortname
                os.makedirs(zsh_startup_dir, exist_ok = True)
                built_zshenv = os.path.join(zsh_startup_dir, '.zshenv')
                shutil.copyfile(os.path.join(srcdir, 'scripts/shells/zsh/zshenv'), built_zshenv)
                built_zshrc = os.path.join(zsh_startup_dir, '.zshrc')
                shutil.copyfile(os.path.join(srcdir, 'scripts/shells/zsh/zshrc.template'), built_zshrc)
                try:
                    fd = open(built_zshrc, 'r+')
                    contents = fd.read()
                    contents = contents.replace('@ENV_VARIABLES@', env_variables)
                    fd.truncate(0)
                    fd.seek(0)
                    fd.write(contents)
                    fd.flush()
                    fd.close()
                except IOError:
                    sys.stderr.write('"{}" failed to update. Check your permissions.\n'.format(built_zshrc))
                    sys.exit(os.EX_CANTCREAT)
        distutils.command.build.build.run(self)

class my_install(install):
    '''
    Override the install to modify updating scripts before installing.
    '''

    def run(self):
        try:
            os.makedirs('build/', exist_ok=True)
        except os.error:
            sys.stderr.write('Build error: failure to create the build/ tree. Please check your permissions.\n')
            sys.exit(os.EX_CANTCREAT)
        # Install is the only time we know the actual data directory.
        # We save this information in a temporary build file for replacement in scripts.
        script = open('build/data_dir', 'w')
        script.truncate(0)
        script.write(os.path.abspath(self.install_data))
        script.close()
        # Go on with normal install.
        install.run(self)

def update_scripts(build_subdir):
    '''
    Convenience function to update any file in `build_subdir`:
    - replace @DATADIR@ by `datadir` as set on the setup.py call.
    '''
    global version
    datadir = '/usr/local'
    # I keep the real version for distribution and release.
    # But if we are in a git repository, the tool will output a git commit too.
    git_version = version
    if shutil.which('git') is not None and os.path.isdir(os.path.join(srcdir, '.git')):
        os.chdir(srcdir)
        if subprocess.check_output(["git", "tag", "--contains"]).decode('utf-8') != 'v' + version:
           commit_hash = subprocess.check_output(['git', 'log', '-1', "--pretty=format:%H"]).decode('utf-8')
           git_version = "development (commit: {} - last release: {})".format(str(commit_hash), version)
        os.chdir(builddir)

    try:
        data_dir_file = open('build/data_dir', 'r')
        datadir = data_dir_file.readline().rstrip(' \n\r\t')
        data_dir_file.close()
    except IOError:
        sys.stderr.write('Warning: no build/data_dir file. You should run the `install` command. Defaulting to {}.\n'.format(datadir))

    for dirpath, dirnames, filenames in os.walk(build_subdir):
        for f in filenames:
            try:
                script = open(os.path.join(dirpath, f), 'r+')
                contents = script.read()
                # Make the necessary replacements.
                contents = contents.replace('@DATADIR@', datadir)
                contents = contents.replace('@VERSION@', git_version)
                script.truncate(0)
                script.seek(0)
                script.write(contents)
                script.flush()
                script.close()
            except IOError:
                sys.stderr.write('The script {} failed to update. Check your permissions.\n'.format(f))
                sys.exit(os.EX_CANTCREAT)

class my_install_data(distutils.command.install_data.install_data):
    '''
    Override the install to build the manual first.
    '''

    def run(self):
        update_scripts('build/platforms')
        update_scripts('build/share/crossroad/scripts')
        distutils.command.install_data.install_data.run(self)
        datadir = '/usr/local'
        try:
            data_dir_file = open('build/data_dir', 'r')
            datadir = data_dir_file.readline().rstrip(' \n\r\t')
            data_dir_file.close()
        except IOError:
            sys.stderr.write('Warning: no build/data_dir file. You should run the `install` command. Defaulting to {}.\n'.format(datadir))
        os.chmod(os.path.join(datadir, 'share/crossroad/scripts/in-crossroad.py'),
                              stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
                              stat.S_IRGRP | stat.S_IXGRP |
                              stat.S_IROTH | stat.S_IXOTH)
        os.chmod(os.path.join(datadir, 'share/crossroad/scripts/config.guess'),
                              stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
                              stat.S_IRGRP | stat.S_IXGRP |
                              stat.S_IROTH | stat.S_IXOTH)
        os.chmod(os.path.join(datadir, 'share/crossroad/scripts/crossroad-mingw-install.py'),
                              stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
                              stat.S_IRGRP | stat.S_IXGRP |
                              stat.S_IROTH | stat.S_IXOTH)
        os.chmod(os.path.join(datadir, 'share/crossroad/scripts/bin-wrappers/crossroad-gcc'),
                              stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
                              stat.S_IRGRP | stat.S_IXGRP |
                              stat.S_IROTH | stat.S_IXOTH)
        os.chmod(os.path.join(datadir, 'share/crossroad/scripts/bin-wrappers/crossroad-cpp'),
                              stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
                              stat.S_IRGRP | stat.S_IXGRP |
                              stat.S_IROTH | stat.S_IXOTH)
        os.chmod(os.path.join(datadir, 'share/crossroad/scripts/bin-wrappers/crossroad-pkg-config'),
                              stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
                              stat.S_IRGRP | stat.S_IXGRP |
                              stat.S_IROTH | stat.S_IXOTH)
        os.makedirs(os.path.join(datadir, 'share/crossroad/bin/'), exist_ok=True)
        # Automatically generate binaries for each platform.
        for f in os.listdir(os.path.join(srcdir, 'platforms/env/')):
            if f[-5:] == '.conf' and f[:-5] not in deactivated_platforms:
                config = configparser.ConfigParser()
                config.read([os.path.join(srcdir, 'platforms/env', f)])
                host = config.get('Platform', 'host')
                if len(host) > 0:
                    try:
                        os.unlink(os.path.join(datadir, 'share/crossroad/bin/' + host + '-pkg-config'))
                    except OSError:
                        pass
                    try:
                        os.unlink(os.path.join(datadir, 'share/crossroad/bin/' + host + '-gcc'))
                    except OSError:
                        pass
                    try:
                        os.unlink(os.path.join(datadir, 'share/crossroad/bin/' + host + '-g++'))
                    except OSError:
                        pass
                    try:
                        os.unlink(os.path.join(datadir, 'share/crossroad/bin/' + host + '-cpp'))
                    except OSError:
                        pass
                    os.symlink(os.path.join(datadir, 'share/crossroad/scripts/bin-wrappers/crossroad-pkg-config'),
                               os.path.join(datadir, 'share/crossroad/bin/' + host + '-pkg-config'))
                    os.symlink(os.path.join(datadir, 'share/crossroad/scripts/bin-wrappers/crossroad-gcc'),
                               os.path.join(datadir, 'share/crossroad/bin/' + host + '-gcc'))
                    os.symlink(os.path.join(datadir, 'share/crossroad/scripts/bin-wrappers/crossroad-gcc'),
                               os.path.join(datadir, 'share/crossroad/bin/' + host + '-g++'))
                    os.symlink(os.path.join(datadir, 'share/crossroad/scripts/bin-wrappers/crossroad-cpp'),
                               os.path.join(datadir, 'share/crossroad/bin/' + host + '-cpp'))

class my_install_scripts(distutils.command.install_scripts.install_scripts):
    '''
    Override the install to build the manual first.
    '''

    def run(self):
        update_scripts(self.build_dir)
        distutils.command.install_scripts.install_scripts.run(self)


platform_list = os.listdir(os.path.join(srcdir, 'platforms/env'))
platform_list = [f[:-5] for f in platform_list \
                 if f[-5:] == '.conf' and f[:-5] not in deactivated_platforms]
platform_file_list = [os.path.join('build/platforms/', f + '.py') for f in platform_list]

cmake_toolchains = [os.path.join('build/share/crossroad/scripts/cmake/',
                                 'toolchain-' + f + '.cmake') \
                    for f in platform_list if f != 'native']
meson_toolchains = [os.path.join('build/share/crossroad/scripts/meson/',
                                 'toolchain-' + f + '.meson') \
                    for f in platform_list if f != 'native']

def get_built_data_files():
    bashrc_files = []
    zsh_files = []
    for f in os.listdir(os.path.join(srcdir, 'platforms/env/')):
        if f[-5:] == '.conf' and f[:-5] not in deactivated_platforms:
            config = configparser.ConfigParser()
            config.read([os.path.join(srcdir, 'platforms/env', f)])
            shortname = config.get('Platform', 'shortname')
            # bash
            built_bashrc = os.path.join('build/share/crossroad/scripts/shells/bash/', 'bashrc.' + shortname)
            bashrc_files += [built_bashrc]
            # ZSH
            zsh_startup_dir = 'build/share/crossroad/scripts/shells/zsh.' + shortname
            built_zshenv = os.path.join(zsh_startup_dir, '.zshenv')
            built_zshrc = os.path.join(zsh_startup_dir, '.zshrc')
            zsh_files += [('share/crossroad/scripts/shells/zsh.' + shortname, [built_zshrc, built_zshenv])]
    return [('share/crossroad/scripts/shells/bash', bashrc_files)] + zsh_files

built_data_files = get_built_data_files()

setup(
    name = 'crossroad',
    cmdclass = {'man': build_man, 'build': my_build, 'install': my_install,
        'install_data': my_install_data, 'install_scripts': my_install_scripts},
    version = version,
    description = 'Cross-Compilation Environment Toolkit.',
    long_description = open(os.path.join(srcdir, 'README')).read(),
    author = 'Jehan',
    author_email = 'jehan@girinstud.io',
    url = 'http://girinstud.io',
    license = 'AGPLv3+',
    classifiers = [
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.3',
        'Topic :: Software Development :: Build Tools',
    ],
    requires = [],
    scripts = ['build/bin/crossroad'],
    data_files = [('share/man/man1/', ['build/share/man/man1/crossroad.1.gz']),
        ('share/crossroad/scripts/', [os.path.join(srcdir, 'scripts/crossroad-mingw-install.py'),
                                      os.path.join(srcdir, 'scripts/config.guess'),
                                      'build/share/crossroad/scripts/in-crossroad.py',
                                      ]),
        ('share/crossroad/scripts/bin-wrappers/',
                                     [os.path.join(srcdir, 'scripts/bin-wrappers/crossroad-gcc'),
                                      os.path.join(srcdir, 'scripts/bin-wrappers/crossroad-pkg-config'),
                                      os.path.join(srcdir, 'scripts/bin-wrappers/crossroad-cpp'),
                                      ]),
        ('share/crossroad/scripts/shells/',
                                     ['build/share/crossroad/scripts/shells/environment.sh',
                                      os.path.join(srcdir, 'scripts/shells/pre-bash-env.sh'),
                                      os.path.join(srcdir, 'scripts/shells/pre-zsh-env.sh'),
                                      os.path.join(srcdir, 'scripts/shells/post-env.sh'),]),
        ('share/bash-completion/completions', [os.path.join(srcdir, 'scripts/shells/bash/completions/crossroad')]),
        ('share/crossroad/scripts/cmake', cmake_toolchains),
        ('share/crossroad/scripts/meson', meson_toolchains),
        ('share/crossroad/platforms/', platform_file_list),
        #('crossroad/projects/', ['projects']),
        ] + built_data_files,
    )