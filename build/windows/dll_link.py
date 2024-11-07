#!/usr/bin/python3

################################################################################
# Small python script to retrieve DLL dependencies with objdump
################################################################################

################################################################################
# Usage example
#
# python3 dll_link.py /path/to/run.exe /winenv/ /path/install
#
# In this case, the DLL dependencies for executable run.exe will be extracted
# and copied into /path/install/bin folder. To copy the DLL, the root path to
# Windows environnment should be passed, here /winenv/.

import argparse
import os
import subprocess
import re
import shutil
import sys

################################################################################
# Global variables

# Sets for executable and system DLL
dlls = set()
sys_dlls = set()

# Common paths
bindir = 'bin'

################################################################################
# Functions

# Main function
def main(binary, srcdir, destdir, debug = False, listfile = None):
    sys.stdout.write("{} (INFO): searching for dependencies of {} in {}\n".format(
        os.path.basename(__file__), binary, srcdir))
    global dlls, sys_dlls
    if listfile is not None:
        listfile = os.path.abspath(listfile)
        if os.path.isfile(listfile):
            print("DLLs listed {:s} will be copied". format(listfile))
            with open(listfile, 'r') as f:
                dll_list = [l.strip() for l in f.readlines()]
                dlls = set(dll_list)
        else:
            print("DLLs list will be saved to {:s}".format(listfile))
            find_dependencies(binary, srcdir)
            installed_dlls = dlls - sys_dlls
            dlldir = os.path.split(listfile)[0]
            if not os.path.isdir(dlldir):
                os.makedirs(dlldir, exist_ok=True)
            with open(listfile, 'w') as f:
                f.writelines(d + "\n" for d in list(installed_dlls))
    else:
        print('Running find_dependencies')
        find_dependencies(binary, srcdir)
    if debug:
        print("Running in debug mode (no DLL moved)")
        if len(dlls) > 0:
            sys.stdout.write("Needed DLLs:\n\t- ")
            sys.stdout.write("\n\t- ".join(list(dlls)))
            print()
        if len(sys_dlls) > 0:
            sys.stdout.write('System DLLs:\n\t- ')
            sys.stdout.write("\n\t- ".join(sys_dlls))
            print()
        removed_dlls = sys_dlls & dlls
        if len(removed_dlls) > 0:
            sys.stdout.write('Non installed DLLs:\n\t- ')
            sys.stdout.write("\n\t- ".join(removed_dlls))
            print()
        installed_dlls = dlls - sys_dlls
        if len(installed_dlls) > 0:
            sys.stdout.write('Installed DLLs:\n\t- ')
            sys.stdout.write("\n\t- ".join(installed_dlls))
            print()
    else:
        copy_dlls(dlls - sys_dlls, srcdir, destdir)


def find_dependencies(obj, srcdir):
    '''
    List DLLs of an object file in a recursive way.
    '''
    if os.path.exists(obj):
        # If DLL exists, extract dependencies.
        objdump = None

        result = subprocess.run(['file', obj], stdout=subprocess.PIPE)
        file_type = result.stdout.decode('utf-8')
        if 'PE32+' in file_type:
            objdump = 'x86_64-w64-mingw32-objdump'
        elif 'PE32' in file_type:
            objdump = 'i686-w64-mingw32-objdump'

        if objdump is None:
            sys.stderr.write(
                'File type of {} unknown: {}\n'.format(obj, file_type))
            sys.exit(os.EX_UNAVAILABLE)
        elif shutil.which(objdump) is None:
            # For native objdump case.
            objdump = 'objdump.exe'

        if shutil.which(objdump) is None:
            sys.stderr.write("Executable doesn't exist: {}\n".format(objdump))
            sys.exit(os.EX_UNAVAILABLE)

        result = subprocess.run([objdump, '-p', obj], stdout=subprocess.PIPE)
        out = result.stdout.decode('utf-8')
        # Parse lines with DLL Name instead of lib*.dll directly
        for match in re.finditer(r"DLL Name: *(\S+.dll)", out, re.MULTILINE):
            dll = match.group(1)
            if dll not in dlls:
                dlls.add(dll)
                next_dll = os.path.join(srcdir, bindir, dll)
                find_dependencies(next_dll, srcdir)
    else:
        # Otherwise, it is a system DLL
        sys_dlls.add(os.path.basename(obj))

# Copy a DLL set into the /destdir/bin directory
def copy_dlls(dll_list, srcdir, destdir):
    destbin = os.path.join(destdir, bindir)
    os.makedirs(destbin, exist_ok=True)
    for dll in dll_list:
        full_file_name = os.path.join(srcdir, bindir, dll)
        if os.path.isfile(full_file_name):
            if not os.path.exists(os.path.join(destbin, dll)):
                sys.stdout.write("{} (INFO): copying {} to {}\n".format(
                    os.path.basename(__file__), full_file_name, destbin))
                shutil.copy(full_file_name, destbin)
        else:
            sys.stderr.write("Missing DLL: {}\n".format(full_file_name))
            sys.exit(os.EX_DATAERR)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', dest='debug',
                        action='store_true', default=False)
    parser.add_argument("--list", "-l", dest = 'listfile', type = str, help="path to list file. If file does not exist it is created, else, it is used to collect the packages", default = None)
    parser.add_argument('binary')
    parser.add_argument('srcdir')
    parser.add_argument('destdir')
    args = parser.parse_args(sys.argv[1:])

    main(**vars(args))
